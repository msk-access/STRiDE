import io
import os

try:
    import torch
    import torch.storage
except ImportError:
    torch = None


def apply_cpu_patches():
    """
    Monkeypatch PyTorch storage reconstructors to force CPU loading of CUDA-saved tensors.
    This avoids the 'Attempting to deserialize object on a CUDA device...' crash on CPU-only runs.
    """
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
