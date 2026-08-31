# Contributing to STRiDE

Thank you for your interest in contributing to STRiDE! This document provides guidelines for contributing to the project.

## Development Setup

```bash
git clone https://github.com/msk-access/STRiDE.git
cd STRiDE
pip install -e ".[dev]"
```

## Git Flow

We follow a Git Flow branching model:

1. Fork or clone the repository
2. Create a feature branch from `develop`: `git checkout -b feature/my-feature develop`
3. Make your changes
4. Run tests: `pytest tests/ -v`
5. Run linting: `ruff check src/ tests/ && black --check src/ tests/`
6. Commit with a descriptive message following [Conventional Commits](https://www.conventionalcommits.org/)
7. Push and open a Pull Request against `develop`

## Code Standards

- **Python**: 3.9+
- **Formatting**: [Black](https://black.readthedocs.io/) (line length 100)
- **Linting**: [Ruff](https://docs.astral.sh/ruff/)
- **Type checking**: [mypy](https://mypy-lang.org/)
- **CLI**: All parameters as `typer.Option()` — no positional arguments
- **Logging**: Use `logging.getLogger(__name__)` — never `print()` or `logging.basicConfig()`
- **Resources**: Access bundled files via `importlib.resources` — never hardcode paths
- **BAM files**: Always use context managers

## Testing

All new features and bug fixes should include tests:

```bash
# Run all tests
pytest tests/ -v

# With coverage
pytest --cov=stride --cov-report=term-missing
```

## Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

- `feat:` — new feature
- `fix:` — bug fix
- `docs:` — documentation changes
- `chore:` — maintenance tasks
- `test:` — adding or updating tests
- `refactor:` — code refactoring without behaviour change

## Reporting Issues

Use [GitHub Issues](https://github.com/msk-access/STRiDE/issues) to report bugs or request features. Include:

- STRiDE version (`stride --version`)
- Python version
- Steps to reproduce
- Expected vs actual behaviour

## License

By contributing, you agree that your contributions will be licensed under the [GNU Affero General Public License v3.0](https://www.gnu.org/licenses/agpl-3.0.html).
