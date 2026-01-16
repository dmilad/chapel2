# Chapel of MOOP - Geodesic Dome Generator

> "Inspired by the principle of Leave No Trace (LNT), this installation promotes playful reflection on our environmental impact through massive art made from MOOP."

**Chapel of MOOP** is a parametric design project for a large-scale art installation proposed for Black Rock City 2026. This repository contains the CadQuery-based design scripts for the structure.

## Main Configuration

The primary dome configuration is defined in `notebooks/chapel.ipynb`:

| Parameter | Value |
|-----------|-------|
| **Style** | Honeycomb (Hex/Pent) |
| **Radius** | 8 ft (243.8 cm) |
| **Frequency** | 3V |
| **Strut Cross-section** | 3.5" × 3.5" |
| **Hub Style** | Tapered Prism |
| **Strut Style** | Cuboid |
| **Windows** | Hexagonal/Pentagonal Plates (2" thick) |

### Generated Part Counts

- **Struts**: 264 total (240 complete, 24 cut at base)
- **Hubs**: 168 total (144 interior, 24 boundary)
- **Windows/Panes**: 89 total (73 complete, 16 partial)

## Project Vision

The **Chapel of MOOP** is a dome-shaped sanctuary constructed with a honeycomb-patterned wooden frame. Its defining feature is the collection of unique windows created by encasing "Matter Out of Place" (MOOP) in dyed epoxy resin.

- **Philosophy**: To draw attention to the LNT principle and our collective environmental impact
- **Dimensions**: 16 feet diameter × 8 feet tall
- **Structure**: Wooden frame with hexagonal/pentagonal cells (honeycomb pattern)
- **Features**:
  - Resin windows encasing collected MOOP
  - LED lighting for nighttime effects
  - Cushioned interior seating

## Project Structure

```
chapel2/
├── notebooks/
│   └── chapel.ipynb      # Main configuration notebook
├── src/chapel2/
│   ├── dome_generator.py # Main dome generation functions
│   ├── geometry.py       # Core geometry calculations
│   ├── analysis.py       # Manufacturability analysis
│   ├── visualization.py  # Selective visualization helpers
│   ├── hubs.py          # Hub joint geometry
│   ├── cuboid_struts.py # Cuboid strut generation
│   ├── wedged_struts.py # Wedged strut generation (alternate style)
│   └── miter_struts.py  # Miter-cut struts (alternate style)
├── output/              # Generated STEP/STL files
├── pyproject.toml       # Poetry configuration
└── README.md
```

## Quick Start

### Prerequisites

- Python 3.12+
- Poetry
- pyenv (recommended)

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

### Usage

Launch JupyterLab and open `notebooks/chapel.ipynb`:

```bash
poetry run jupyter lab
```

For visualization, use VS Code/Cursor with the OCP CAD Viewer extension:
1. Open VS Code/Cursor
2. Press `Cmd+Shift+P` → "OCP CAD Viewer: Open Viewer"
3. Run notebook cells with `show()` to visualize geometry

## Key Functions

### Dome Generation

```python
from chapel2.dome_generator import generate_dome_with_hubs

struts, hubs, windows, info = generate_dome_with_hubs(
    radius_cm=8.0 * 30.48,      # 8 ft in cm
    frequency=3,
    strut_width=3.5 * 2.54,     # 3.5" in cm
    strut_depth=3.5 * 2.54,
    dome_style="honeycomb",
    hub_style="tapered_prism",
    strut_style="cuboid",
    generate_windows=True,
    window_plate_depth=2 * 2.54, # 2" thick
)
```

### Selective Visualization

```python
from chapel2.visualization import separate_dome_parts, VISUALIZATION_COLORS

separated = separate_dome_parts(
    vertices=info['vertices'],
    edges=info['edges'],
    faces=info['faces'],
    radius_cm=info['radius_cm'],
    portion=info['portion'],
    strut_width=info['strut_width'],
    strut_depth=info['strut_depth'],
    hub_inset=info['hub_inset'],
    window_plate_depth=2 * 2.54,
    window_margin=0.2,
)

# Returns dict with:
# - complete_struts, partial_struts
# - complete_hubs, boundary_hubs
# - complete_panes, partial_panes
```

### Manufacturability Analysis

```python
from chapel2.analysis import full_manufacturability_analysis, format_analysis_report

analysis = full_manufacturability_analysis(
    vertices, edges, faces,
    radius_cm, portion, strut_width, hub_inset
)
print(format_analysis_report(analysis))
```

## Configuration Options

### Hub Styles
- `tapered_prism` - Flat triangular inner face (default)
- `convex_hull` - Wraps strut corners in convex hull
- `cylindrical_core` - Cylinder for miter-cut struts

### Strut Styles
- `cuboid` - Rectangular cross-section (default)
- `wedged` - Trapezoidal with radial side faces
- `miter_cut` - For cylindrical core hubs

### Dome Styles
- `honeycomb` - Hex/pent pattern (default)
- `triangular` - Standard triangular geodesic

## Credits

- **Lead Artist**: Milad
- **Event**: Burning Man 2026 (Honoraria Letter of Intent)

## Additional Notebooks

The `notebooks/` directory also contains development notebooks:
- `dome_test.ipynb` - Testing various dome configurations
- `hub_comparison.ipynb` - Comparing hub styles
- `example.ipynb` - Basic CadQuery examples

These are kept for reference but the main configuration is in `chapel.ipynb`.
