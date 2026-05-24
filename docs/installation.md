# Installation

## Without installing: run directly with uv

If you just want to run a script that uses `eagar-tsai` without adding it to a project or creating a virtual environment, uv can handle that in one step:

```bash
uv run --with eagar-tsai python my_script.py
```

uv fetches the package into a temporary, isolated environment and runs your script inside it. Nothing is installed globally or into any project. This is the fastest way to try the package or run a one-off calculation on a machine where you do not want to manage a project environment.

## With uv (add to a project)

```bash
uv add eagar-tsai
```

## With pip

```bash
pip install eagar-tsai
```

Pre-built binary wheels are published to PyPI for Python 3.11, 3.12, and 3.13 on Linux (x86-64, i686), macOS (x86-64 and Apple Silicon), and Windows (AMD64). If a matching wheel exists for your platform, no C compiler is needed; pip or uv will download and install the pre-compiled package directly.

!!! warning "Building from source"
    If no pre-built wheel matches your platform (for example, a non-standard Linux architecture or a Python version outside the supported range), the package falls back to building from the source distribution. In that case a C compiler is required:

    - **Linux/macOS**: GCC or Clang (usually pre-installed)
    - **Windows**: MSVC Build Tools or MinGW-w64
