"""
CadQuery Geodesic Dome Generator with Wedged Struts.

This module generates geodesic dome structures with wedged struts - trapezoidal 
prisms whose side faces are radial planes passing through the dome center.
"""

import cadquery as cq
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace
from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCP.gp import gp_Pnt
from OCP.TopoDS import TopoDS_Shape
from typing import List, Tuple, Dict, Optional

from .geometry import (
    Point3D,
    Edge,
    generate_geodesic_dome,
    generate_honeycomb_dome,
    extract_polygon_edges,
    normalize,
    cross,
    add,
    sub,
    scale,
    norm,
)


# =============================================================================
# WEDGED STRUT GEOMETRY
# =============================================================================

def project_radially(point: Point3D, distance: float, center: Point3D) -> Point3D:
    """
    Project a point radially outward from the dome center by a given distance.
    
    Args:
        point: The point to project
        distance: How far to project outward (positive = away from center)
        center: The dome center point
    
    Returns:
        The projected point
    """
    radial = normalize(sub(point, center))
    return add(point, scale(radial, distance))


def compute_profile_corners(
    p: Point3D, 
    other: Point3D, 
    width: float, 
    depth: float, 
    dome_center: Point3D
) -> List[Point3D]:
    """
    Compute the 4 corners of a strut profile at a given endpoint.
    
    The strut cross-section is a trapezoid whose side faces are radial planes
    (they pass through the dome center).
    
    Args:
        p: The endpoint where we're computing the profile
        other: The other endpoint of the strut (defines the strut axis)
        width: Width of the strut (tangential direction)
        depth: Depth of the strut (radial direction, inner to outer)
        dome_center: Center point of the dome
    
    Returns:
        List of 4 corners: [inner_right, inner_left, outer_left, outer_right]
        in a clockwise order when viewed from outside the dome
    """
    # Radial direction at point p (pointing outward from dome center)
    radial = normalize(sub(p, dome_center))
    
    # Strut axis direction (from p toward other)
    axis = sub(other, p)
    
    # Side normal (perpendicular to both radial and axis)
    # This gives us the "width" direction of the strut
    side_normal = normalize(cross(radial, axis))
    
    # Handle degenerate case where radial and axis are parallel
    if norm(side_normal) < 1e-9:
        # Fallback to an arbitrary perpendicular
        if abs(radial[2]) < 0.9:
            side_normal = normalize(cross(radial, (0, 0, 1)))
        else:
            side_normal = normalize(cross(radial, (1, 0, 0)))
    
    half_w = width / 2.0
    
    # Inner corners (at the current radius, on the inner surface)
    inner_left = add(p, scale(side_normal, half_w))
    inner_right = sub(p, scale(side_normal, half_w))
    
    # Outer corners (projected radially outward by depth)
    outer_left = project_radially(inner_left, depth, dome_center)
    outer_right = project_radially(inner_right, depth, dome_center)
    
    return [inner_right, inner_left, outer_left, outer_right]


def create_wedged_strut(
    corners_start: List[Point3D], 
    corners_end: List[Point3D]
) -> cq.Shape:
    """
    Create a wedged strut solid using CadQuery/OCP's loft (ThruSections).
    
    Args:
        corners_start: 4 corner points at the start end
        corners_end: 4 corner points at the end (must be in matching order)
    
    Returns:
        CadQuery Shape representing the strut solid
    """
    # Use OCP's ThruSections to loft between two wire profiles
    builder = BRepOffsetAPI_ThruSections(True, False)  # True = solid, False = ruled
    
    # Create start wire (closed polygon)
    wire_start = BRepBuilderAPI_MakePolygon()
    for pt in corners_start:
        wire_start.Add(gp_Pnt(pt[0], pt[1], pt[2]))
    wire_start.Close()
    builder.AddWire(wire_start.Wire())
    
    # Create end wire (closed polygon)
    wire_end = BRepBuilderAPI_MakePolygon()
    for pt in corners_end:
        wire_end.Add(gp_Pnt(pt[0], pt[1], pt[2]))
    wire_end.Close()
    builder.AddWire(wire_end.Wire())
    
    builder.Build()
    
    return cq.Shape(builder.Shape())


def create_strut_from_edge(
    start: Point3D,
    end: Point3D,
    strut_width: float,
    strut_depth: float,
    dome_center: Point3D = (0, 0, 0)
) -> cq.Shape:
    """
    Create a single wedged strut between two vertices.
    
    Args:
        start: Start vertex position
        end: End vertex position
        strut_width: Width of the strut cross-section
        strut_depth: Depth of the strut (radial thickness)
        dome_center: Center of the dome
    
    Returns:
        CadQuery Shape of the strut
    """
    # Compute profile corners at each end
    corners_start = compute_profile_corners(start, end, strut_width, strut_depth, dome_center)
    corners_end = compute_profile_corners(end, start, strut_width, strut_depth, dome_center)
    
    # Reorder end corners to match start corners for proper loft
    # The end profile is computed looking back toward start, so we need to flip it
    corners_end = [corners_end[1], corners_end[0], corners_end[3], corners_end[2]]
    
    return create_wedged_strut(corners_start, corners_end)


# =============================================================================
# DOME GENERATION
# =============================================================================

def generate_dome(
    radius_cm: float,
    frequency: int,
    portion: float = 0.5,
    strut_width: float = 5.08,
    strut_depth: float = 5.08,
    dome_style: str = "honeycomb",
    max_struts: int = -1
) -> Tuple[cq.Compound, Dict]:
    """
    Generate a complete geodesic dome with wedged struts.
    
    Args:
        radius_cm: Dome radius in centimeters
        frequency: Geodesic frequency (1-6, higher = more subdivisions)
        portion: Sphere portion (0.5 = hemisphere)
        strut_width: Width of strut cross-section in cm
        strut_depth: Depth of strut (radial thickness) in cm
        dome_style: "honeycomb" (hex/pent) or "triangular"
        max_struts: Maximum struts to generate (-1 = all)
    
    Returns:
        Tuple of (compound solid, info dict with vertices/edges/faces)
    """
    dome_center = (0.0, 0.0, 0.0)
    
    # Generate dome geometry based on style
    if dome_style == "honeycomb":
        vertices, faces, edges = generate_honeycomb_dome(radius_cm, frequency, portion)
    else:
        vertices, faces, edges = generate_geodesic_dome(radius_cm, frequency, portion)
    
    # Limit struts if requested
    if max_struts > 0:
        edges = edges[:max_struts]
    
    # Create all struts
    strut_shapes = []
    
    for v1_idx, v2_idx in edges:
        start = vertices[v1_idx]
        end = vertices[v2_idx]
        
        strut = create_strut_from_edge(start, end, strut_width, strut_depth, dome_center)
        strut_shapes.append(strut)
    
    # Combine all struts into one compound
    if strut_shapes:
        result = cq.Compound.makeCompound(strut_shapes)
    else:
        result = cq.Compound.makeCompound([])
    
    info = {
        'num_vertices': len(vertices),
        'num_edges': len(edges),
        'num_faces': len(faces),
        'vertices': vertices,
        'edges': edges,
        'faces': faces,
        'radius_cm': radius_cm,
        'frequency': frequency,
        'portion': portion,
        'strut_width': strut_width,
        'strut_depth': strut_depth,
    }
    
    return result, info


def generate_dome_struts_individually(
    radius_cm: float,
    frequency: int,
    portion: float = 0.5,
    strut_width: float = 5.08,
    strut_depth: float = 5.08,
    dome_style: str = "honeycomb",
    max_struts: int = -1
) -> Tuple[List[cq.Shape], Dict]:
    """
    Generate dome struts as individual shapes (useful for debugging/visualization).
    
    Same parameters as generate_dome, but returns a list of individual strut shapes
    instead of a compound.
    """
    dome_center = (0.0, 0.0, 0.0)
    
    # Generate dome geometry based on style
    if dome_style == "honeycomb":
        vertices, faces, edges = generate_honeycomb_dome(radius_cm, frequency, portion)
    else:
        vertices, faces, edges = generate_geodesic_dome(radius_cm, frequency, portion)
    
    # Limit struts if requested
    if max_struts > 0:
        edges = edges[:max_struts]
    
    # Create all struts
    strut_shapes = []
    
    for v1_idx, v2_idx in edges:
        start = vertices[v1_idx]
        end = vertices[v2_idx]
        
        strut = create_strut_from_edge(start, end, strut_width, strut_depth, dome_center)
        strut_shapes.append(strut)
    
    info = {
        'num_vertices': len(vertices),
        'num_edges': len(edges),
        'num_faces': len(faces),
        'vertices': vertices,
        'edges': edges,
        'faces': faces,
        'radius_cm': radius_cm,
        'frequency': frequency,
        'portion': portion,
        'strut_width': strut_width,
        'strut_depth': strut_depth,
    }
    
    return strut_shapes, info


# =============================================================================
# EXPORT FUNCTIONS
# =============================================================================

def export_dome(
    dome: cq.Compound,
    step_path: Optional[str] = None,
    stl_path: Optional[str] = None
) -> None:
    """
    Export the dome to STEP and/or STL files.
    
    Args:
        dome: The dome compound to export
        step_path: Path for STEP file (None to skip)
        stl_path: Path for STL file (None to skip)
    """
    # Wrap the compound in a Workplane for export
    wp = cq.Workplane("XY").newObject([dome])
    
    if step_path:
        cq.exporters.export(wp, step_path)
        print(f"Exported STEP to: {step_path}")
    
    if stl_path:
        cq.exporters.export(wp, stl_path)
        print(f"Exported STL to: {stl_path}")


# =============================================================================
# CONVENIENCE / PRESET FUNCTIONS
# =============================================================================

def generate_2v_test_dome(
    radius_cm: float = 100.0,
    strut_width: float = 5.0,
    strut_depth: float = 5.0,
    max_struts: int = 10
) -> Tuple[cq.Compound, Dict]:
    """
    Generate a simple 2V dome for testing.
    
    This is useful for verifying the geometry before generating a full dome.
    """
    return generate_dome(
        radius_cm=radius_cm,
        frequency=2,
        portion=0.5,
        strut_width=strut_width,
        strut_depth=strut_depth,
        dome_style="honeycomb",
        max_struts=max_struts
    )


def generate_chapel_dome(
    radius_ft: float = 8.0,
    frequency: int = 3,
    strut_width_in: float = 2.0,
    strut_depth_in: float = 2.0
) -> Tuple[cq.Compound, Dict]:
    """
    Generate the Chapel of MOOP dome with standard parameters.
    
    Args:
        radius_ft: Dome radius in feet
        frequency: Geodesic frequency (3 recommended)
        strut_width_in: Strut width in inches
        strut_depth_in: Strut depth in inches
    
    Returns:
        Tuple of (compound, info dict)
    """
    # Convert to centimeters
    radius_cm = radius_ft * 30.48
    strut_width = strut_width_in * 2.54
    strut_depth = strut_depth_in * 2.54
    
    return generate_dome(
        radius_cm=radius_cm,
        frequency=frequency,
        portion=0.5,
        strut_width=strut_width,
        strut_depth=strut_depth,
        dome_style="honeycomb"
    )


# =============================================================================
# MAIN (for standalone testing)
# =============================================================================

if __name__ == "__main__":
    print("Generating test dome...")
    
    # Generate a simple 2V test dome
    dome, info = generate_2v_test_dome(radius_cm=100, max_struts=10)
    
    print(f"Generated dome with:")
    print(f"  Vertices: {info['num_vertices']}")
    print(f"  Edges: {info['num_edges']}")
    print(f"  Faces: {info['num_faces']}")
    
    # Export to files
    export_dome(dome, step_path="output/test_dome.step", stl_path="output/test_dome.stl")
    
    print("Done!")

