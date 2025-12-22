# Design Reference

This document captures the key configuration options and parameters from the original Fusion 360 implementation.

## Dome Structure Types

The chapel supports two dome styles:
- **Triangular**: Classic geodesic dome (Class I) - triangular faces
- **Honeycomb**: Dual mesh with hexagons and pentagons - **The Chapel of MOOP design**

## Preset Configurations

### 1. 2ft Hexagons (3V Dome)
- **Radius**: 6.74 ft (205.4 cm)
- **Frequency**: 3V
- **Hexagons**: 2.00 ft side-to-side
- **Pentagons**: 1.61 ft side-to-side
- **Structure**: 89 faces, 286 struts

### 2. 8ft Dome (3V)
- **Radius**: 8.00 ft (243.84 cm)
- **Frequency**: 3V
- **Hexagons**: 2.37 ft side-to-side
- **Pentagons**: 1.91 ft side-to-side
- **Structure**: 89 faces, 286 struts

### 3. 4V 8ft Dome (Higher Detail)
- **Radius**: 8.00 ft (243.84 cm)
- **Frequency**: 4V
- **Hexagons**: 1.20 ft side-to-side
- **Pentagons**: 0.96 ft side-to-side
- **Structure**: 337 faces, 1054 struts

## Main Parameters

### Dome Dimensions
- **DOME_RADIUS**: Radius in centimeters (default: 205.4 cm / 6.74 ft)
  - Common conversions:
    - 8 feet = 243.84 cm
    - 10 feet = 304.8 cm
    - 12 feet = 365.76 cm
    - 16 feet = 487.68 cm

- **FREQUENCY**: Geodesic subdivision level
  - 1V = 20 triangles (basic icosahedron)
  - 2V = 80 triangles
  - 3V = 180 triangles
  - 4V = 320 triangles
  - Higher frequency = smoother dome but more struts to build

- **DOME_PORTION**: Fraction of sphere to include
  - 0.5 = hemisphere (flat bottom)
  - 0.6 = slightly taller than hemisphere
  - 1.0 = full sphere

### Strut Dimensions (Picture Frame Style)
- **STRUT_WIDTH**: Cuboid width (default: 5.08 cm / 2")
- **STRUT_DEPTH**: Cuboid depth/thickness (default: 5.08 cm / 2")
- Common nominal lumber sizes:
  - 2" × 2" = 5.08 cm
  - 3" × 3" = 7.62 cm
  - 4" × 4" = 10.16 cm

**Note**: Real lumber is often smaller than nominal (e.g., "2×2" ≈ 1.5"×1.5")

### Hub Connectors
- **CREATE_SOLID_HUBS**: Enable solid triangular hub connectors (default: True)
- **HUB_INSET**: How much to shorten each strut end for hub joints (default: STRUT_WIDTH / 2.0)
  - Recommended: ~half of STRUT_WIDTH for clean joints

## Current Target: 16ft Diameter × 8ft Tall

For a 16ft diameter dome:
- **Radius**: 8 ft = 243.84 cm
- **Height**: 8 ft (requires calculating DOME_PORTION)
- **Recommended Frequency**: 3V or 4V
  - 3V: 89 faces, ~2.37 ft hexagons, 286 struts
  - 4V: 337 faces, ~1.20 ft hexagons, 1054 struts

To achieve 8ft height with 8ft radius:
- For a hemisphere (DOME_PORTION = 0.5), height = radius = 8ft ✓

## Geometry Generation Process

1. **Create Icosahedron**: Start with 12 vertices, 20 triangular faces
2. **Subdivide**: Apply frequency-based subdivision (each triangle → 4 triangles per level)
3. **Normalize**: Project vertices onto unit sphere
4. **Compute Dual** (for honeycomb): Convert triangular mesh to hexagon/pentagon mesh
5. **Scale**: Apply radius
6. **Filter**: Keep only upper portion (dome)
7. **Extract Edges**: Generate strut list

## Window Analysis

The honeycomb dome produces:
- **Hexagons**: 6-sided polygons (majority of windows)
- **Pentagons**: 5-sided polygons (12 total in a full sphere, ~6-10 in hemisphere)

Window dimensions include:
- **Outer side-to-side**: Frame outer dimension
- **Inner side-to-side**: Actual window opening (accounting for strut width)
- **Edge lengths**: For calculating resin panel sizes

## Assembly Outputs

The design system should generate:
1. **Cut List**: Grouped strut lengths with counts
2. **Hub Types**: Different joint angles require different hub geometries
3. **Window Dimensions**: Sizes for resin panels (hexagons and pentagons)
4. **Assembly Guide**: Vertex-to-hub mapping, strut connections
5. **3D Model**: For visualization and fabrication planning
6. **Export Formats**: STEP, STL, DXF for manufacturing
