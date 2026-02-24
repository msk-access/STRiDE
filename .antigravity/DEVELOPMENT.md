# STRiDE Development Guide

## Environment Setup

```bash
# Clone and install in editable mode with dev dependencies
git clone https://github.com/msk-access/STRiDE.git
cd STRiDE
pip install -e ".[dev]"

# Verify installation
stride --version
pytest tests/ -v
```

## Git Flow Branching

```
main          ← stable releases (tagged)
  └── develop      ← integration branch (default)
       └── feature/*   ← individual features/fixes
```

- All work happens on `feature/*` branches created from `develop`
- PRs merge into `develop`; `develop` merges into `main` for releases
- Never push directly to `main`

## Release Process

1. Create release branch: `git checkout -b release/X.Y.Z develop`
2. Bump `version` in `pyproject.toml` and `src/stride/__init__.py`
3. Update `CHANGELOG.md` (move Unreleased → version header)
4. Commit, push, PR to `main`
5. Merge to `main`, tag: `git tag X.Y.Z && git push origin X.Y.Z`
6. CI automatically: builds → publishes to PyPI → pushes Docker to GHCR
7. Merge tag back to `develop`

## Testing

```bash
# Run all tests
pytest tests/ -v

# With coverage
pytest --cov=stride --cov-report=term-missing

# CLI-specific tests use typer.testing.CliRunner
```

## Code Conventions

- **CLI**: All parameters as `typer.Option()` — no positional arguments
- **Logging**: Use `logging.getLogger(__name__)` in every module; never call `logging.basicConfig()` in library code
- **Rich**: `setup_logging()` in `utils.py` is the single logging entry point
- **Resources**: Access bundled files via `importlib.resources.files()`, never hardcode paths
- **BAMs**: Always use context managers (`with MSIProfileGenerator(...) as gen:`)

## CI/CD Workflows

| File | Trigger | Purpose |
|------|---------|---------|
| `test.yml` | Push/PR → develop/main | Test matrix + lint + Docker build |
| `release.yml` | Tag push | PyPI + GHCR publish |
| `deploy-docs.yml` | Push (docs paths) | MkDocs versioned deploy via mike |
