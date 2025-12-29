# Chapel of MOOP (CadQuery Edition)

> "Inspired by the principle of Leave No Trace (LNT), this installation promotes playful reflection on our environmental impact through massive art made from MOOP."

**Chapel of MOOP** is a parametric design project for a large-scale art installation proposed for Black Rock City 2026. This repository contains the CadQuery-based design scripts for the structure.

## Project Vision

The **Chapel of MOOP** is a dome-shaped sanctuary constructed with a honeycomb-patterned wooden frame. Its defining feature is the collection of unique windows created by encasing "Matter Out of Place" (MOOP) in dyed epoxy resin.

- **Philosophy**: To draw attention to the LNT principle and our collective environmental impact. Inside the chapel, surrounded by MOOP-filled stained glass, participants are invited to reflect on how small, discarded objects accumulate into a larger footprint.
- **Physical Description**:
  - **Dimensions**: 16 feet diameter × 8 feet tall (configurable).
  - **Structure**: Wooden frame with hexagonal/pentagonal cells (honeycomb pattern).
  - **Frame Style**: Picture-frame struts (rectangular cross-section) with hub connectors at joints. Supports multiple hub styles:
    - Tapered triangular prism hubs (flat inner face, tapered outer)
    - Cylindrical core hubs with miter-cut struts
    - Convex hull hubs (baseline)
  - **Features**:
    - Windows: Dyed resin encasing collected MOOP (with some protruding for tactile interaction).
    - Lighting: LED lights surrounding windows for nighttime patterns.
    - Interior: Cushions for seating, hidden chest with MOOP-made gifts.

## The Design Tool

This project uses CadQuery (via Jupyter notebooks) to parametrically generate the chapel's geometry. This allows for rapid iteration on:
- Dome frequency (3V, 4V) and radius
- Hexagon/pentagon window sizing (e.g., 2ft side-to-side)
- Picture-frame strut dimensions (width × depth)
- Solid triangular hub connectors at joints
- Automatic cut list generation

The scripts generate full 3D geometry including rectangular struts and hub connectors, ready for fabrication planning and export to STEP/STL formats.

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
├── notebooks/          # Jupyter notebooks for dome design and modeling
│   └── example.ipynb   # CadQuery examples
├── src/
│   └── chapel2/       # Python modules (geometry, parameters, etc.)
├── pyproject.toml     # Poetry configuration
└── README.md
```

## Credits

- **Lead Artist**: Milad
- **Event**: Burning Man 2026 (Honoraria Letter of Intent)

## Notes

This is a reimplementation of the Chapel of MOOP design using CadQuery instead of Fusion 360. The original Fusion 360-based implementation can be found in the [chapel](https://github.com/dmilad/chapel) repository.
