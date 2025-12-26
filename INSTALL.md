# Installation Guide

## Prerequisites

Before you can use `bonded_substructures`, you need to install the following system dependencies.

## Step 1: Install Python Package Management Tools

### Option A: Using apt (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install -y python3-pip python3-venv
```

### Option B: Install pip manually (if sudo not available)

Download and run the pip installer:

```bash
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py --user
```

Add pip to your PATH:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

## Step 2: Install Poetry (Recommended)

Poetry provides better dependency management and virtual environment handling:

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

Add Poetry to your PATH:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

Verify installation:

```bash
poetry --version
```

## Step 3: Install Project Dependencies

### Option A: Using Poetry (Recommended)

```bash
cd bonded_substructures
poetry install
```

This creates a virtual environment and installs all dependencies.

### Option B: Using pip with venv

```bash
cd bonded_substructures
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e .
```

### Option C: Install dependencies manually

If you just want to install the core dependencies:

```bash
pip install --user gmsh numpy scipy matplotlib pyvista
pip install --user fenics-dolfinx
```

For development:

```bash
pip install --user pytest pytest-cov black ruff
```

## Step 4: Verify Installation

### Using Poetry

```bash
poetry run python examples/01_rectangle_basic.py
```

### Using pip/venv

```bash
python examples/01_rectangle_basic.py
```

## System-Specific Dependencies

### gmsh

On some systems, you may need to install gmsh system-wide:

**Ubuntu/Debian:**
```bash
sudo apt-get install gmsh
```

**macOS:**
```bash
brew install gmsh
```

### dolfinx (FEniCSx)

dolfinx has additional system dependencies. See the [FEniCSx installation guide](https://docs.fenicsproject.org/dolfinx/main/python/installation.html).

**Ubuntu/Debian (recommended):**
```bash
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenicsx
```

**Using conda (alternative):**
```bash
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx
```

## Troubleshooting

### "No module named 'gmsh'"

Install gmsh:
```bash
pip install --user gmsh
```

### "No module named 'dolfinx'"

Follow the FEniCSx installation guide above. Note that dolfinx is only needed for future FE analysis features, not for basic mesh generation.

### "poetry: command not found"

Make sure Poetry is in your PATH:
```bash
export PATH="$HOME/.local/bin:$PATH"
```

### Permission denied errors

Use `--user` flag with pip:
```bash
pip install --user <package>
```

Or use a virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate
pip install <package>
```

## Quick Start After Installation

Once installed, try running the examples:

```bash
# Using Poetry
poetry run python examples/01_rectangle_basic.py
poetry run python examples/02_rectangle_disbond.py

# Or activate the Poetry shell
poetry shell
python examples/01_rectangle_basic.py

# Using venv
source venv/bin/activate
python examples/01_rectangle_basic.py
```

## Running Tests

```bash
# Using Poetry
poetry run pytest

# Using venv
source venv/bin/activate
pytest
```
