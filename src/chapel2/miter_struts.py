"""
Miter-Cut Strut Geometry for Geodesic Domes.

Miter-cut struts have compound angle cuts at the ends so the end faces 
lie in radial planes (planes that pass through the dome center).
This allows struts to meet at a simple cylindrical hub core.
"""

import math
import cadquery as cq
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCP.gp import gp_Pnt
from typing import List, Tuple, Dict, Optional

from .geometry import (
    Point3D,
    normalize,
    sub,
    add,
    scale,
    cross,
    dot,
    norm
)


def compute_miter_cut_angles(
    strut_direction: Point3D,
    radial_at_vertex: Point3D
) -> Dict[str, float]:
    """
    Compute the compound miter angles for cutting a strut end.
    
    The goal is to cut the strut end so the face lies in a "radial plane" -
    a plane that contains the radial direction at the vertex.
    
    For woodworking/manufacturing:
    - Miter angle: rotation around the strut's depth axis (in the tangent plane)
    - Bevel angle: rotation around the strut's width axis (tilting in/out)
    
    Args:
        strut_direction: Unit vector along the strut (toward the vertex)
        radial_at_vertex: Unit vector pointing radially outward at the vertex
    
    Returns:
        Dict with:
        - 'miter_angle_deg': Miter angle in degrees
        - 'bevel_angle_deg': Bevel angle in degrees  
        - 'cut_plane_normal': Normal to the cut plane
    """
    # The strut comes in at some angle to the radial direction
    # A normal perpendicular cut would have the face perpendicular to strut_direction
    # We want to rotate this face to lie in the radial plane
    
    # The radial plane at the vertex contains:
    # 1. The radial direction
    # 2. The strut direction projected onto the tangent plane
    
    # Compute the angle between strut direction and radial
    cos_angle = dot(strut_direction, radial_at_vertex)
    strut_to_radial_angle = math.acos(max(-1, min(1, cos_angle)))
    
    # The cut plane normal should be in the plane spanned by 
    # strut_direction and radial_at_vertex
    # Specifically, it should be perpendicular to the line where
    # the strut meets the radial plane
    
    # Side normal (perpendicular to both radial and strut direction)
    side_normal = normalize(cross(radial_at_vertex, strut_direction))
    
    # Handle degenerate case where strut is exactly radial
    if norm(side_normal) < 1e-9:
        return {
            'miter_angle_deg': 0.0,
            'bevel_angle_deg': 0.0,
            'cut_plane_normal': strut_direction,
        }
    
    # The cut plane should contain:
    # - The side_normal (so the width of the strut lies in the cut face)
    # - A direction in the radial plane
    
    # For a proper miter cut, we bisect the angle between the strut
    # and the "outgoing" direction to the vertex center
    # But for radial plane cuts, we want the cut face to contain the radial direction
    
    # The cut normal is perpendicular to both:
    # 1. side_normal (already in tangent plane)
    # 2. radial direction
    # This means cut_normal = normalize(side_normal × radial)
    # But this would be parallel to the tangent plane projection of strut direction!
    
    # Actually, we want the cut face to CONTAIN the radial direction
    # So the cut plane normal is perpendicular to radial
    # And also perpendicular to side_normal
    # cut_normal = radial × side_normal
    
    cut_normal = normalize(cross(radial_at_vertex, side_normal))
    
    # The miter angle is the angle between this cut_normal and the strut_direction
    # (both projected onto the tangent plane if needed)
    miter_cos = dot(cut_normal, strut_direction)
    miter_angle_rad = math.acos(max(-1, min(1, abs(miter_cos))))
    
    # The bevel angle relates to the tilt out of the tangent plane
    # For cuboid struts with tangential width, bevel is typically 0
    bevel_angle_rad = 0.0
    
    return {
        'miter_angle_deg': math.degrees(miter_angle_rad),
        'bevel_angle_deg': math.degrees(bevel_angle_rad),
        'cut_plane_normal': cut_normal,
        'strut_to_radial_angle_deg': math.degrees(strut_to_radial_angle),
    }


def compute_miter_cut_corners(
    vertex: Point3D,
    other: Point3D,
    strut_width: float,
    strut_depth: float,
    hub_inset: float,
    dome_center: Point3D
) -> List[Point3D]:
    """
    Compute the 4 corners of a miter-cut strut end face.
    
    The cut face lies in a radial plane (contains the radial direction),
    which allows it to meet a cylindrical hub tangentially.
    
    For a proper miter cut:
    - The inner edge (toward dome center) is CLOSER to the vertex
    - The outer edge is FARTHER from the vertex
    - This creates an angled cut face that lies in a radial plane
    
    Args:
        vertex: The vertex where this strut end meets the hub
        other: The other end of the strut
        strut_width: Width of the strut
        strut_depth: Depth of the strut (radial thickness)
        hub_inset: How far the strut CENTER is shortened from the vertex
        dome_center: Center of the dome
    
    Returns:
        List of 4 corners of the miter-cut end face
    """
    # Radial direction at vertex (pointing outward from dome center)
    radial = normalize(sub(vertex, dome_center))
    
    # Strut direction (from other toward vertex - direction strut is coming from)
    strut_dir = normalize(sub(vertex, other))
    
    # Side normal (tangential, defines width direction)
    # This is perpendicular to both radial and strut direction
    side_normal = normalize(cross(radial, strut_dir))
    if norm(side_normal) < 1e-9:
        # Strut is exactly radial - use arbitrary perpendicular
        if abs(radial[2]) < 0.9:
            side_normal = normalize(cross(radial, (0, 0, 1)))
        else:
            side_normal = normalize(cross(radial, (1, 0, 0)))
    
    half_w = strut_width / 2.0
    half_d = strut_depth / 2.0
    
    # Compute the angle between strut direction and radial
    # cos(theta) = strut_dir . radial
    cos_theta = dot(strut_dir, radial)
    
    # For the cut face to lie in a radial plane:
    # - The face normal must be perpendicular to radial
    # - The face must contain points at different distances along the strut
    #   depending on their radial position (inner vs outer)
    
    # The key insight: if we move by 'half_d' in the radial direction,
    # we need to move by 'half_d / tan(theta)' along the strut direction
    # to stay in the same radial plane.
    #
    # But we need sin(theta) for the perpendicular component
    sin_theta = math.sqrt(max(0, 1 - cos_theta * cos_theta))
    
    # Avoid division by zero for nearly-radial struts
    if sin_theta < 0.1:
        # Strut is nearly radial - use small angle approximation
        # The miter cut is nearly perpendicular
        strut_offset_per_radial = 0.0
    else:
        # How much to move along strut axis per unit of radial offset
        # to keep the point in the same radial plane
        strut_offset_per_radial = cos_theta / sin_theta
    
    # Base point: strut end center at hub_inset from vertex
    # This is the CENTER of the cut face
    end_center = add(vertex, scale(strut_dir, -hub_inset))
    
    # For inner corners (toward dome center, offset by -half_d in radial):
    # Need to move CLOSER to vertex (less negative along strut_dir)
    inner_strut_offset = -hub_inset + half_d * strut_offset_per_radial
    inner_center = add(vertex, scale(strut_dir, inner_strut_offset))
    inner_center = sub(inner_center, scale(radial, half_d))
    
    # For outer corners (away from dome center, offset by +half_d in radial):
    # Need to move FARTHER from vertex (more negative along strut_dir)
    outer_strut_offset = -hub_inset - half_d * strut_offset_per_radial
    outer_center = add(vertex, scale(strut_dir, outer_strut_offset))
    outer_center = add(outer_center, scale(radial, half_d))
    
    # Now add width offsets along side_normal
    corners = [
        sub(inner_center, scale(side_normal, half_w)),  # inner right
        add(inner_center, scale(side_normal, half_w)),  # inner left
        add(outer_center, scale(side_normal, half_w)),  # outer left
        sub(outer_center, scale(side_normal, half_w)),  # outer right
    ]
    
    return corners


def create_miter_cut_strut(
    start: Point3D,
    end: Point3D,
    strut_width: float,
    strut_depth: float,
    hub_inset_start: float,
    hub_inset_end: float,
    dome_center: Point3D = (0, 0, 0)
) -> cq.Shape:
    """
    Create a cuboid strut with miter-cut ends.
    
    Both ends are cut so the end faces lie in radial planes,
    allowing them to meet cylindrical hub cores.
    
    Args:
        start: Start vertex position
        end: End vertex position
        strut_width: Width of the strut
        strut_depth: Depth of the strut
        hub_inset_start: How far to shorten at start end
        hub_inset_end: How far to shorten at end
        dome_center: Center of the dome
    
    Returns:
        CadQuery Shape of the strut
    """
    # Get corners for both ends
    corners_start = compute_miter_cut_corners(
        start, end, strut_width, strut_depth, hub_inset_start, dome_center
    )
    corners_end = compute_miter_cut_corners(
        end, start, strut_width, strut_depth, hub_inset_end, dome_center
    )
    
    # Reorder end corners to match start for proper loft
    # The corners_end are computed from end's perspective, need to flip
    corners_end = [corners_end[1], corners_end[0], corners_end[3], corners_end[2]]
    
    # Loft between the two end profiles
    builder = BRepOffsetAPI_ThruSections(True, False)  # True = solid, False = ruled
    
    # Create start wire
    wire_start = BRepBuilderAPI_MakePolygon()
    for pt in corners_start:
        wire_start.Add(gp_Pnt(pt[0], pt[1], pt[2]))
    wire_start.Close()
    builder.AddWire(wire_start.Wire())
    
    # Create end wire
    wire_end = BRepBuilderAPI_MakePolygon()
    for pt in corners_end:
        wire_end.Add(gp_Pnt(pt[0], pt[1], pt[2]))
    wire_end.Close()
    builder.AddWire(wire_end.Wire())
    
    builder.Build()
    
    return cq.Shape(builder.Shape())


def analyze_miter_cuts_for_dome(
    vertices: List[Point3D],
    edges: List[Tuple[int, int]],
    dome_center: Point3D = (0, 0, 0),
    angle_tolerance_deg: float = 0.5
) -> Dict:
    """
    Analyze all miter cuts needed for a dome and group by similar angles.
    
    This helps determine how many unique strut types need to be manufactured.
    
    Args:
        vertices: All vertex positions
        edges: All edges (vertex index pairs)
        dome_center: Center of the dome
        angle_tolerance_deg: How close angles must be to group together
    
    Returns:
        Dict with cut analysis including unique cut configurations
    """
    all_cuts = []
    
    for edge_idx, (v1_idx, v2_idx) in enumerate(edges):
        v1 = vertices[v1_idx]
        v2 = vertices[v2_idx]
        
        # Compute miter angles at both ends
        radial_v1 = normalize(sub(v1, dome_center))
        radial_v2 = normalize(sub(v2, dome_center))
        
        strut_dir_from_v1 = normalize(sub(v2, v1))
        strut_dir_from_v2 = normalize(sub(v1, v2))
        
        cut_v1 = compute_miter_cut_angles(strut_dir_from_v1, radial_v1)
        cut_v2 = compute_miter_cut_angles(strut_dir_from_v2, radial_v2)
        
        # Strut length
        length = norm(sub(v2, v1))
        
        all_cuts.append({
            'edge_idx': edge_idx,
            'length': length,
            'cut_angle_start_deg': cut_v1['miter_angle_deg'],
            'cut_angle_end_deg': cut_v2['miter_angle_deg'],
        })
    
    # Group struts by similar (length, cut_angle_start, cut_angle_end) tuples
    def round_to_tolerance(val, tol):
        return round(val / tol) * tol
    
    strut_groups = {}
    for cut in all_cuts:
        key = (
            round_to_tolerance(cut['length'], 0.1),  # 0.1cm = 1mm tolerance
            round_to_tolerance(cut['cut_angle_start_deg'], angle_tolerance_deg),
            round_to_tolerance(cut['cut_angle_end_deg'], angle_tolerance_deg),
        )
        if key not in strut_groups:
            strut_groups[key] = []
        strut_groups[key].append(cut)
    
    return {
        'total_struts': len(edges),
        'unique_strut_types': len(strut_groups),
        'strut_groups': strut_groups,
        'all_cuts': all_cuts,
    }


