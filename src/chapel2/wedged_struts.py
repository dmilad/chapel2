"""
Wedged Strut Geometry for Geodesic Domes.

Wedged struts are trapezoidal prisms whose side faces are radial planes 
passing through the dome center.
"""

import cadquery as cq
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCP.gp import gp_Pnt
from typing import List

from .geometry import (
    Point3D,
    normalize,
    sub,
    add,
    scale,
    cross,
    norm,
    project_radially
)


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





