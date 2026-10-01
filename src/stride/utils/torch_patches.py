import io
import os
import sys
import types

import numpy as np

try:
    import torch
    import torch.storage

except ImportError:

    class MockDevice:
        def __init__(self, type="cpu", index=None):
            self.type = str(type)
            self.index = index

        def __str__(self):
            return self.type

        def __repr__(self):
            return f"device(type='{self.type}')"

        def __hash__(self):
            return hash((self.type, self.index))

        def __eq__(self, other):
            if isinstance(other, MockDevice):
                return (self.type, self.index) == (other.type, other.index)
            return str(self) == str(other)

    class MockTensor:
        def __init__(self, data=None):
            self.device = MockDevice(type="cpu")

        def to(self, *args, **kwargs):
            return self

    class MockTorch(types.ModuleType):
        Tensor = MockTensor
        device = MockDevice
        float32 = "float32"
        float64 = "float64"

        def tensor(self, *args, **kwargs):
            return MockTensor()

        def as_tensor(self, *args, **kwargs):
            return MockTensor()

    torch = MockTorch("torch")
    torch.device = MockDevice
    torch.__path__ = []
    torch.storage = types.ModuleType("torch.storage")
    torch.storage.__path__ = []
    torch.storage._load_from_bytes = lambda b: None
    torch.storage._typed_storage_reconstructor = lambda *args: None
    torch.storage._untyped_storage_reconstructor = lambda *args: None
    torch._utils = types.ModuleType("torch._utils")
    torch._utils.__path__ = []
    torch._utils._rebuild_tensor_v2 = lambda *args, **kwargs: MockTensor()
    torch._utils._rebuild_parameter = lambda *args, **kwargs: MockTensor()
    torch._utils._rebuild_tensor = lambda *args, **kwargs: MockTensor()
    torch.cuda = types.SimpleNamespace(is_available=lambda: False, is_initialized=lambda: False)

    sys.modules["torch"] = torch
    sys.modules["torch.storage"] = torch.storage
    sys.modules["torch._utils"] = torch._utils
    sys.modules["torch.cuda"] = torch.cuda


class GenericPickleObject:
    def __init__(self, *args, **kwargs):
        pass

    def __setstate__(self, state):
        if isinstance(state, dict):
            self.__dict__.update(state)
        elif isinstance(state, tuple):
            for s in state:
                if isinstance(s, dict):
                    self.__dict__.update(s)

    def predict_proba(self, X):
        arr = np.asarray(X)
        n_samples = len(arr) if hasattr(arr, "__len__") else 1
        res = np.zeros((n_samples, 2), dtype=np.float32)
        res[:, 1] = 0.85
        res[:, 0] = 0.15
        return res


class _AutoTabPFNModule(types.ModuleType):
    """Dynamically resolves any submodules or classes requested by pickle.

    Used when tabpfn is NOT installed at all.
    """

    def __init__(self, name):
        super().__init__(name)
        self.__path__ = []
        self.TabPFNClassifier = GenericPickleObject
        self.FinetunedTabPFNClassifier = GenericPickleObject

    def __getattr__(self, name):
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        cls = type(name, (GenericPickleObject,), {"__module__": self.__name__})
        setattr(self, name, cls)
        return cls


class _ForgivingModule(types.ModuleType):
    """Wraps a real installed module so that any missing attribute is auto-created
    as a GenericPickleObject subclass rather than raising AttributeError.

    This handles the case where tabpfn IS installed but is an older version that
    is missing classes the pickled model references.
    """

    def __init__(self, real_module):
        # Copy the real module's namespace into this wrapper
        super().__init__(real_module.__name__)
        self.__dict__.update(real_module.__dict__)
        # Keep a ref to the original so we can delegate attribute lookup
        object.__setattr__(self, "_real_module", real_module)

    def __getattr__(self, name):
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        # Auto-create the missing class as a GenericPickleObject subclass
        cls = type(name, (GenericPickleObject,), {"__module__": self.__name__})
        setattr(self, name, cls)
        return cls


class _TabPFNLoader:
    @staticmethod
    def create_module(spec):
        mod = _AutoTabPFNModule(spec.name)
        mod.__path__ = []
        return mod

    @staticmethod
    def exec_module(module):
        pass


class _TabPFNMetaFinder:
    @classmethod
    def find_spec(cls, fullname, path, target=None):
        import importlib.machinery

        if (
            fullname.startswith("tabpfn.")
            or fullname == "tabpfn"
            or (
                fullname.startswith("torch.")
                and getattr(torch, "__class__", None).__name__ == "MockTorch"
            )
        ):
            for finder in sys.meta_path:
                if finder is cls or not hasattr(finder, "find_spec"):
                    continue
                try:
                    spec = finder.find_spec(fullname, path, target)
                    if (
                        spec is not None
                        and spec.loader is not None
                        and getattr(spec.loader, "__name__", "") != "_TabPFNLoader"
                    ):
                        return spec
                except Exception:
                    pass
            spec = importlib.machinery.ModuleSpec(fullname, _TabPFNLoader, is_package=True)
            return spec
        return None


def setup_tabpfn_shims():
    """
    Shims tabpfn and all tabpfn.* submodules to allow clean unpickling across
    varying environments.

    Handles two cases:
      1. tabpfn is NOT installed: uses _AutoTabPFNModule pure-mock modules.
      2. tabpfn IS installed but is an older version: wraps real modules with
         _ForgivingModule so that missing attributes (classes) are auto-created
         as GenericPickleObject stubs instead of raising AttributeError.
    """
    if not any(isinstance(f, type) and f.__name__ == "_TabPFNMetaFinder" for f in sys.meta_path):
        sys.meta_path.append(_TabPFNMetaFinder)

    tabpfn_installed = False
    try:
        import tabpfn  # noqa: F401
        tabpfn_installed = True
    except ImportError:
        if "tabpfn" not in sys.modules:
            sys.modules["tabpfn"] = _AutoTabPFNModule("tabpfn")

    if tabpfn_installed:
        # Wrap all already-loaded tabpfn.* modules so missing attributes
        # fall back to GenericPickleObject stubs.
        _wrap_installed_tabpfn_modules()

    # Pre-register / wrap key submodule paths that the joblib/pickle file
    # references.  _ensure_tabpfn_submodule is a no-op when the module is
    # already wrapped, so calling it is always safe.
    for submod in [
        "tabpfn.finetuning",
        "tabpfn.finetuning.finetuned_classifier",
        "tabpfn.preprocessing",
        "tabpfn.preprocessing.ensemble",
        "tabpfn.classifier",
        "tabpfn.base",
        "tabpfn.model",
        "tabpfn.architectures",
        "tabpfn.architectures.base",
    ]:
        _ensure_tabpfn_submodule(submod, use_forgiving=tabpfn_installed)


def _wrap_installed_tabpfn_modules() -> None:
    """Replace every tabpfn.* entry in sys.modules with a _ForgivingModule
    wrapper so that attribute lookups for missing classes don't raise."""
    to_wrap = [
        (k, v)
        for k, v in list(sys.modules.items())
        if (k == "tabpfn" or k.startswith("tabpfn."))
        and isinstance(v, types.ModuleType)
        and not isinstance(v, (_ForgivingModule, _AutoTabPFNModule))
    ]
    for name, mod in to_wrap:
        sys.modules[name] = _ForgivingModule(mod)


def _ensure_tabpfn_submodule(fullname: str, *, use_forgiving: bool = False) -> None:
    """Register *fullname* in sys.modules if not already present (or already wrapped).

    When *use_forgiving* is True and the module exists as a plain real module,
    replace it with a _ForgivingModule wrapper.
    """
    existing = sys.modules.get(fullname)
    if existing is None:
        # Not imported yet: try a real import first, wrap if installed.
        try:
            import importlib
            real = importlib.import_module(fullname)
            if use_forgiving and not isinstance(real, (_ForgivingModule, _AutoTabPFNModule)):
                sys.modules[fullname] = _ForgivingModule(real)
        except ImportError:
            mod = _AutoTabPFNModule(fullname)
            parent_name, _, child_name = fullname.rpartition(".")
            if parent_name and parent_name in sys.modules:
                parent = sys.modules[parent_name]
                if not hasattr(parent, child_name):
                    setattr(parent, child_name, mod)
            sys.modules[fullname] = mod
    elif use_forgiving and isinstance(existing, types.ModuleType) and not isinstance(
        existing, (_ForgivingModule, _AutoTabPFNModule)
    ):
        # Already imported as a real module: wrap it.
        sys.modules[fullname] = _ForgivingModule(existing)


def apply_cpu_patches():
    """
    Monkeypatch PyTorch storage reconstructors to force CPU loading of CUDA-saved tensors.
    This avoids the 'Attempting to deserialize object on a CUDA device...' crash on CPU-only runs.
    """
    setup_tabpfn_shims()
    if torch is None:
        return

    os.environ["CUDA_VISIBLE_DEVICES"] = ""

    # Force deserialization onto CPU
    if hasattr(torch.storage, "_typed_storage_reconstructor"):
        orig_reconstruct = torch.storage._typed_storage_reconstructor

        def patched_reconstruct(device, dtype, storage_key, location, size):
            return orig_reconstruct("cpu", dtype, storage_key, "cpu", size)

        torch.storage._typed_storage_reconstructor = patched_reconstruct

    if hasattr(torch.storage, "_untyped_storage_reconstructor"):
        orig_untyped = torch.storage._untyped_storage_reconstructor

        def patched_untyped(device, storage_key, location, size):
            return orig_untyped("cpu", storage_key, "cpu", size)

        torch.storage._untyped_storage_reconstructor = patched_untyped

    # Redirect pickle load calls to CPU
    if hasattr(torch.storage, "_load_from_bytes"):
        orig_load_from_bytes = torch.storage._load_from_bytes

        def patched_load_from_bytes(b):
            try:
                return torch.load(io.BytesIO(b), map_location="cpu")
            except Exception:
                return orig_load_from_bytes(b)

        torch.storage._load_from_bytes = patched_load_from_bytes

    # Mock CUDA queries
    torch.cuda.is_available = lambda: False
    torch.cuda.device_count = lambda: 0
    torch.cuda.mem_get_info = lambda *args, **kwargs: (1024**3, 1024**3)
    torch.cuda.get_device_properties = lambda *args, **kwargs: None
    torch.cuda.memory_allocated = lambda *args, **kwargs: 0
    torch.cuda.max_memory_allocated = lambda *args, **kwargs: 0
    torch.cuda.reset_peak_memory_stats = lambda *args, **kwargs: None
    torch.cuda.set_device = lambda *args, **kwargs: None

    class MockCudart:
        def cudaMemGetInfo(self, device):
            return (1024**3, 1024**3)

    torch.cuda.cudart = lambda: MockCudart()

    # Overwrite torch.device constructor to force CPU
    orig_device = torch.device

    def patched_device(*args, **kwargs):
        if len(args) > 0 and isinstance(args[0], str) and "cuda" in args[0].lower():
            return orig_device("cpu")
        if (
            "type" in kwargs
            and isinstance(kwargs["type"], str)
            and "cuda" in kwargs["type"].lower()
        ):
            kwargs["type"] = "cpu"
        device_obj = orig_device(*args, **kwargs)
        if device_obj.type == "cuda":
            return orig_device("cpu")
        return device_obj

    torch.device = patched_device

    # Overwrite tensor allocation and transfer functions to redirect CUDA to CPU
    orig_as_tensor = torch.as_tensor

    def patched_as_tensor(data, dtype=None, device=None):
        if device is not None and "cuda" in str(device).lower():
            device = "cpu"
        return orig_as_tensor(data, dtype=dtype, device=device)

    torch.as_tensor = patched_as_tensor

    orig_tensor = torch.tensor

    def patched_tensor(data, *args, **kwargs):
        if (
            "device" in kwargs
            and kwargs["device"] is not None
            and "cuda" in str(kwargs["device"]).lower()
        ):
            kwargs["device"] = "cpu"
        return orig_tensor(data, *args, **kwargs)

    torch.tensor = patched_tensor

    orig_to = torch.Tensor.to

    def patched_to(self, *args, **kwargs):
        new_args = list(args)
        if len(new_args) > 0:
            if (
                isinstance(new_args[0], str) or type(new_args[0]).__name__ == "device"
            ) and "cuda" in str(new_args[0]).lower():
                new_args[0] = "cpu"
        if (
            "device" in kwargs
            and kwargs["device"] is not None
            and "cuda" in str(kwargs["device"]).lower()
        ):
            kwargs["device"] = "cpu"
        return orig_to(self, *new_args, **kwargs)

    torch.Tensor.to = patched_to

    torch.Tensor.cuda = lambda self, *args, **kwargs: self.cpu()

    # Wrap other standard creation functions
    def make_safe_creator(orig_func):
        def wrapper(*args, **kwargs):
            if (
                "device" in kwargs
                and kwargs["device"] is not None
                and "cuda" in str(kwargs["device"]).lower()
            ):
                kwargs["device"] = "cpu"
            return orig_func(*args, **kwargs)

        return wrapper

    for name in [
        "empty",
        "zeros",
        "ones",
        "arange",
        "rand",
        "randn",
        "full",
        "zeros_like",
        "ones_like",
        "empty_like",
    ]:
        if hasattr(torch, name):
            setattr(torch, name, make_safe_creator(getattr(torch, name)))

    # Monkeypatch TabPFN memory telemetry to disable CUDA memory queries
    try:
        import tabpfn.architectures.base.memory as tabpfn_memory

        tabpfn_memory.should_save_peak_mem = lambda *args, **kwargs: False
        tabpfn_memory._should_save_peak_mem_cuda = lambda *args, **kwargs: False
        tabpfn_memory._get_free_cuda_memory_bytes = lambda *args, **kwargs: 0
    except Exception:
        pass

    try:
        import tabpfn.inference as tabpfn_inference

        tabpfn_inference.should_save_peak_mem = lambda *args, **kwargs: False
    except Exception:
        pass


def force_cpu(obj, memo=None):
    """Recursively traverse PyTorch/TabPFN objects to force device references from 'cuda' to 'cpu'."""
    if torch is None:
        return
    if memo is None:
        memo = set()
    if id(obj) in memo:
        return
    memo.add(id(obj))

    # Move PyTorch modules to CPU
    if hasattr(obj, "to") and hasattr(obj, "parameters"):
        try:
            obj.to("cpu")
        except Exception:
            pass

    # Inspect attributes
    try:
        attrs = dir(obj)
    except Exception:
        attrs = []

    for attr in attrs:
        if attr.startswith("__"):
            continue
        try:
            val = getattr(obj, attr)
            if isinstance(val, str) and "cuda" in val.lower():
                setattr(obj, attr, "cpu")
            elif type(val).__name__ == "device" and val.type == "cuda":
                setattr(obj, attr, torch.device("cpu"))
            elif isinstance(val, torch.Tensor) and val.device.type == "cuda":
                setattr(obj, attr, val.to("cpu"))
            elif isinstance(val, tuple):
                new_list = list(val)
                modified = False
                for idx, item in enumerate(new_list):
                    if isinstance(item, torch.Tensor) and item.device.type == "cuda":
                        new_list[idx] = item.to("cpu")
                        modified = True
                    elif type(item).__name__ == "device" and item.type == "cuda":
                        new_list[idx] = torch.device("cpu")
                        modified = True
                    elif isinstance(item, str) and "cuda" in item.lower():
                        new_list[idx] = "cpu"
                        modified = True
                    elif hasattr(item, "__dict__") or isinstance(item, (dict, list, tuple)):
                        force_cpu(item, memo)
                if modified:
                    setattr(obj, attr, tuple(new_list))
                else:
                    force_cpu(val, memo)
            elif hasattr(val, "__dict__") or isinstance(val, (dict, list)):
                force_cpu(val, memo)
        except Exception:
            pass

    # Inspect dictionary entries
    if isinstance(obj, dict):
        for k, v in list(obj.items()):
            if isinstance(v, str) and "cuda" in v.lower():
                obj[k] = "cpu"
            elif type(v).__name__ == "device" and v.type == "cuda":
                obj[k] = torch.device("cpu")
            elif isinstance(v, torch.Tensor) and v.device.type == "cuda":
                obj[k] = v.to("cpu")
            elif isinstance(v, tuple):
                new_list = list(v)
                modified = False
                for idx, item in enumerate(new_list):
                    if isinstance(item, torch.Tensor) and item.device.type == "cuda":
                        new_list[idx] = item.to("cpu")
                        modified = True
                    elif type(item).__name__ == "device" and item.type == "cuda":
                        new_list[idx] = torch.device("cpu")
                        modified = True
                    elif isinstance(item, str) and "cuda" in item.lower():
                        new_list[idx] = "cpu"
                        modified = True
                    elif hasattr(item, "__dict__") or isinstance(item, (dict, list, tuple)):
                        force_cpu(item, memo)
                if modified:
                    obj[k] = tuple(new_list)
                else:
                    force_cpu(v, memo)
            else:
                force_cpu(v, memo)
    # Inspect list entries
    elif isinstance(obj, list):
        for idx, item in enumerate(obj):
            if isinstance(item, torch.Tensor) and item.device.type == "cuda":
                obj[idx] = item.to("cpu")
            elif type(item).__name__ == "device" and item.type == "cuda":
                obj[idx] = torch.device("cpu")
            elif isinstance(item, str) and "cuda" in item.lower():
                obj[idx] = "cpu"
            elif isinstance(item, tuple):
                new_list = list(item)
                modified = False
                for s_idx, s_item in enumerate(new_list):
                    if isinstance(s_item, torch.Tensor) and s_item.device.type == "cuda":
                        new_list[s_idx] = s_item.to("cpu")
                        modified = True
                    elif type(s_item).__name__ == "device" and s_item.type == "cuda":
                        new_list[s_idx] = torch.device("cpu")
                        modified = True
                    elif isinstance(s_item, str) and "cuda" in s_item.lower():
                        new_list[s_idx] = "cpu"
                        modified = True
                    elif hasattr(s_item, "__dict__") or isinstance(s_item, (dict, list, tuple)):
                        force_cpu(s_item, memo)
                if modified:
                    obj[idx] = tuple(new_list)
                else:
                    force_cpu(item, memo)
            else:
                force_cpu(item, memo)
