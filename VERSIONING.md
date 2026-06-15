# Versioning and releases

This repository uses [`bump-my-version`](https://callowayproject.github.io/bump-my-version/) to keep release identifiers synchronized across the Python package, documentation, citation metadata, and project summary.

Current release: **0.1.0**

## Install release tooling

```bash
python -m pip install -e ".[dev]"
```

The `dev` extra installs `bump-my-version` together with the test dependencies.

## Inspect the current version

```bash
bump-my-version show current_version
bump-my-version show-bump
```

## Create a version bump

Run these commands from a clean Git working tree:

```bash
bump-my-version bump patch
bump-my-version bump minor
bump-my-version bump major
```

The configuration in `pyproject.toml` updates:

- `pyproject.toml`
- `src/ga2o3dichroism/__init__.py`
- `README.md`
- `CITATION.cff`
- `NOTICE`
- `VERSIONING.md`
- `package-summary.json`

The release bump also creates a commit with the message:

```text
chore(release): bump version <old> → <new>
```

and a Git tag named:

```text
v<new>
```

## Recommended release check

```bash
pytest -q
python -m compileall src tests
bump-my-version show current_version
```
