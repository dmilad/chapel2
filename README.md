# Chapel2

CadQuery project for 3D modeling with Jupyter visualization.

## Setup

This project uses Poetry for dependency management and pyenv for Python version management.

### Prerequisites

- Python 3.12+
- Poetry
- pyenv

### Installation

1. Clone the repository
2. Set up the virtual environment:
   ```bash
   pyenv virtualenv 3.12.12 chapel2
   pyenv local chapel2
   ```

3. Install dependencies:
   ```bash
   poetry install
   ```

## Usage

Launch JupyterLab to start creating 3D models:

```bash
poetry run jupyter lab
```

Open the example notebook in the `notebooks/` directory to see CadQuery in action.

## Project Structure

```
chapel2/
├── notebooks/          # Jupyter notebooks for modeling
├── src/
│   └── chapel2/       # Python modules
├── pyproject.toml     # Poetry configuration
└── README.md
```
