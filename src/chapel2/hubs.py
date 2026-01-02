"""
Hub Joint Geometry for Geodesic Domes.

Hubs are connectors at vertices where struts meet.
"""

import math
import cadquery as cq
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakePolygon, 
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeSolid,
    BRepBuilderAPI_Sewing
)
from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCP.gp import gp_Pnt
from OCP.TopoDS import TopoDS, TopoDS_Shape
from typing import List, Dict, Optional, Tuple

from .geometry import (
    Point3D,
    Edge,
    normalize,
    sub,
    add,
    scale,
    cross,
    norm,
    dot,
    compute_intersection_2d,
    project_to_tangent_plane,
    compute_tangent_basis,
    point_3d_to_2d,
    point_2d_to_3d
)


def calculate_hub_corner(
    vertex: Point3D,
    radial: Point3D,
    d1: Point3D,
    d2: Point3D,
    hub_inset: float,
    strut_depth: float
) -> Point3D:
    """
    Calculate the corner of a hub between two struts.
    This corresponds to the intersection of the two strut cut planes on the inner face.
    """
    # Tangent basis at vertex
    u, v = compute_tangent_basis(radial)
    
    # Hub inner face center (inset by half depth)
    inner_center = sub(vertex, scale(radial, strut_depth / 2.0))
    
    # For each strut, finding the cut plane line on the inner face
    # The cut plane passes through (vertex + hub_inset * direction)
    # And is perpendicular to direction
    
    # We work in 2D on the inner face plane
    
    # Strut 1
    p1_3d = add(vertex, scale(d1, hub_inset))
    p1_2d = point_3d_to_2d(p1_3d, inner_center, u, v)
    
    # Direction 1 projected
    d1_proj = project_to_tangent_plane(d1, radial)
    d1_2d = (
        d1_proj[0]*u[0] + d1_proj[1]*u[1] + d1_proj[2]*u[2],
        d1_proj[0]*v[0] + d1_proj[1]*v[1] + d1_proj[2]*v[2]
    )
    
    # Line 1 direction is perpendicular to d1_2d
    l1_dir = (-d1_2d[1], d1_2d[0])
    
    # Strut 2
    p2_3d = add(vertex, scale(d2, hub_inset))
    p2_2d = point_3d_to_2d(p2_3d, inner_center, u, v)
    
    d2_proj = project_to_tangent_plane(d2, radial)
    d2_2d = (
        d2_proj[0]*u[0] + d2_proj[1]*u[1] + d2_proj[2]*u[2],
        d2_proj[0]*v[0] + d2_proj[1]*v[1] + d2_proj[2]*v[2]
    )
    
    l2_dir = (-d2_2d[1], d2_2d[0])
    
    # Intersection
    corner_2d = compute_intersection_2d(p1_2d, l1_dir, p2_2d, l2_dir)
    
    if corner_2d is None:
        # Parallel lines? Should not happen for adjacent struts
        return inner_center
        
    return point_2d_to_3d(corner_2d, inner_center, u, v)


def compute_hub_geometry(
    vertex: Point3D,
    vertex_idx: int,
    vertices: List[Point3D],
    edges: List[Edge],
    vertex_to_edges: Dict[int, List[int]],
    strut_width: float,
    strut_depth: float,
    hub_inset: float,
    dome_center: Point3D = (0, 0, 0)
) -> Dict:
    """
    Compute the geometry for a hub at a given vertex.
    
    Args:
        vertex: The vertex position
        vertex_idx: Index of this vertex
        vertices: All vertices
        edges: All edges
        vertex_to_edges: Map from vertex index to edge indices
        strut_width: Width of struts
        strut_depth: Depth of struts (radial)
        hub_inset: How far struts are shortened from vertex
        dome_center: Center of dome
    
    Returns:
        Dict with hub geometry info:
        - 'vertex': vertex position
        - 'strut_directions': list of unit vectors pointing toward connected vertices
        - 'strut_endpoints': list of points where shortened struts would end
        - 'num_struts': number of struts meeting at this vertex
    """
    edge_indices = vertex_to_edges.get(vertex_idx, [])
    
    strut_directions = []
    strut_endpoints = []
    connected_vertices = []
    
    for edge_idx in edge_indices:
        v1_idx, v2_idx = edges[edge_idx]
        # Find the other vertex
        other_idx = v2_idx if v1_idx == vertex_idx else v1_idx
        other_vertex = vertices[other_idx]
        
        # Direction from this vertex toward the other
        direction = normalize(sub(other_vertex, vertex))
        strut_directions.append(direction)
        connected_vertices.append(other_idx)
        
        # Endpoint where shortened strut would end (inset from vertex)
        endpoint = add(vertex, scale(direction, hub_inset))
        strut_endpoints.append(endpoint)
    
    return {
        'vertex': vertex,
        'vertex_idx': vertex_idx,
        'strut_directions': strut_directions,
        'strut_endpoints': strut_endpoints,
        'connected_vertices': connected_vertices,
        'num_struts': len(edge_indices),
        'hub_inset': hub_inset,
    }


def create_hub_solid(
    hub_info: Dict,
    strut_width: float,
    strut_depth: float,
    dome_center: Point3D = (0, 0, 0)
) -> Optional[cq.Shape]:
    """
    Create a triangular prism hub at a vertex where struts meet.
    
    For 3-strut vertices: creates a triangular prism with faces perpendicular 
    to each incoming strut direction.
    
    For 2-strut vertices: creates a simple bar connector.
    
    Args:
        hub_info: Hub geometry from compute_hub_geometry()
        strut_width: Width of struts
        strut_depth: Depth of struts
        dome_center: Center of dome
    
    Returns:
        CadQuery Shape for the hub, or None if invalid
    """
    vertex = hub_info['vertex']
    strut_directions = hub_info['strut_directions']
    hub_inset = hub_info['hub_inset']
    num_struts = hub_info['num_struts']
    
    if num_struts < 2:
        return None
    
    # Radial direction at the vertex (this is the prism axis direction)
    radial = normalize(sub(vertex, dome_center))
    
    if num_struts == 3:
        # Create convex hull for 3-strut vertices (default)
        return _create_convex_hull_hub_3strut(
            vertex, radial, strut_directions, hub_inset, strut_width, strut_depth
        )
    elif num_struts == 2:
        # Create a simple bar connector for 2-strut (edge) vertices
        return _create_bar_hub(
            vertex, radial, strut_directions, hub_inset, strut_width, strut_depth
        )
    else:
        # For vertices with more than 3 struts, use convex hull approach
        return _create_convex_hull_hub(
            hub_info, strut_width, strut_depth, dome_center
        )


def _compute_edge_intersection_2d(
    p1: Tuple[float, float], d1: Tuple[float, float],
    p2: Tuple[float, float], d2: Tuple[float, float]
) -> Optional[Tuple[float, float]]:
    """Legacy wrapper for compute_intersection_2d."""
    return compute_intersection_2d(p1, d1, p2, d2)


def _project_to_tangent_plane(
    direction: Point3D,
    radial: Point3D
) -> Point3D:
    """Legacy wrapper for project_to_tangent_plane."""
    return project_to_tangent_plane(direction, radial)


def _compute_tangent_basis(radial: Point3D) -> Tuple[Point3D, Point3D]:
    """Legacy wrapper for compute_tangent_basis."""
    return compute_tangent_basis(radial)


def _point_3d_to_2d(point: Point3D, origin: Point3D, u: Point3D, v: Point3D) -> Tuple[float, float]:
    """Legacy wrapper for point_3d_to_2d."""
    return point_3d_to_2d(point, origin, u, v)


def _point_2d_to_3d(point_2d: Tuple[float, float], origin: Point3D, u: Point3D, v: Point3D) -> Point3D:
    """Legacy wrapper for point_2d_to_3d."""
    return point_2d_to_3d(point_2d, origin, u, v)


def _create_tapered_prism_hub(
    vertex: Point3D,
    radial: Point3D,
    strut_directions: List[Point3D],
    hub_inset: float,
    strut_width: float,
    strut_depth: float
) -> Optional[cq.Shape]:
    """
    Create a tapered triangular prism hub for 3-strut vertices.
    
    The inner face is a flat triangle in the tangent plane (perpendicular to radial).
    The outer face is another triangle, offset along the radial direction.
    The three side faces connect corresponding edges of the triangles.
    
    This design provides a clean flat inner surface for aesthetics.
    """
    if len(strut_directions) != 3:
        # Fall back to convex hull for non-3-strut vertices
        return None
    
    half_w = strut_width / 2.0
    half_d = strut_depth / 2.0
    
    # Set up tangent plane coordinate system
    u, v = _compute_tangent_basis(radial)
    
    # Struts are centered on the geodesic edge (vertex position)
    # Half extends outward, half extends inward. Hub matches this.
    outer_center = add(vertex, scale(radial, half_d))  # Outer face at half_d above vertex
    inner_center = sub(vertex, scale(radial, half_d))  # Inner face at half_d below vertex
    
    # For each strut, compute the inner edge line in 2D tangent plane coordinates
    # The inner edge is at distance hub_inset along strut direction, offset by half_width perpendicular
    
    inner_edge_lines = []  # List of (point_2d, direction_2d) for inner edges
    outer_edge_lines = []  # List of (point_2d, direction_2d) for outer edges
    
    for direction in strut_directions:
        # Project strut direction onto tangent plane
        d_tangent = _project_to_tangent_plane(direction, radial)
        
        # Side normal in tangent plane (perpendicular to projected direction)
        side_normal = normalize(cross(radial, direction))
        side_tangent = _project_to_tangent_plane(side_normal, radial)
        
        # The inner edge of this strut at the hub:
        # Start at vertex, move hub_inset along direction, then this is the strut end face
        # The inner edge is at the inner face of the strut (toward dome center)
        
        # Point on the inner edge (center of the inner edge of strut end face)
        inner_edge_center_3d = add(vertex, scale(direction, hub_inset))
        inner_edge_center_2d = _point_3d_to_2d(inner_edge_center_3d, inner_center, u, v)
        
        # Direction along the inner edge (perpendicular to strut direction in tangent plane)
        edge_dir_2d = (
            side_tangent[0]*u[0] + side_tangent[1]*u[1] + side_tangent[2]*u[2],
            side_tangent[0]*v[0] + side_tangent[1]*v[1] + side_tangent[2]*v[2]
        )
        
        inner_edge_lines.append((inner_edge_center_2d, edge_dir_2d))
        
        # For outer edges, same logic but offset outward
        outer_edge_center_2d = _point_3d_to_2d(inner_edge_center_3d, outer_center, u, v)
        outer_edge_lines.append((outer_edge_center_2d, edge_dir_2d))
    
    # Compute inner triangle vertices (intersection of adjacent edge lines)
    inner_triangle_2d = []
    outer_triangle_2d = []
    
    for i in range(3):
        j = (i + 1) % 3
        
        # Inner triangle vertex = intersection of edge i and edge j
        inner_pt = _compute_edge_intersection_2d(
            inner_edge_lines[i][0], inner_edge_lines[i][1],
            inner_edge_lines[j][0], inner_edge_lines[j][1]
        )
        if inner_pt is None:
            return None  # Degenerate case
        inner_triangle_2d.append(inner_pt)
        
        # Outer triangle vertex
        outer_pt = _compute_edge_intersection_2d(
            outer_edge_lines[i][0], outer_edge_lines[i][1],
            outer_edge_lines[j][0], outer_edge_lines[j][1]
        )
        if outer_pt is None:
            return None
        outer_triangle_2d.append(outer_pt)
    
    # Convert back to 3D
    inner_triangle_3d = [_point_2d_to_3d(pt, inner_center, u, v) for pt in inner_triangle_2d]
    outer_triangle_3d = [_point_2d_to_3d(pt, outer_center, u, v) for pt in outer_triangle_2d]
    
    # Build the tapered prism by lofting between inner and outer triangles
    try:
        builder = BRepOffsetAPI_ThruSections(True, True)  # Solid, ruled
        
        # Inner triangle wire
        wire_inner = BRepBuilderAPI_MakePolygon()
        for pt in inner_triangle_3d:
            wire_inner.Add(gp_Pnt(pt[0], pt[1], pt[2]))
        wire_inner.Close()
        builder.AddWire(wire_inner.Wire())
        
        # Outer triangle wire
        wire_outer = BRepBuilderAPI_MakePolygon()
        for pt in outer_triangle_3d:
            wire_outer.Add(gp_Pnt(pt[0], pt[1], pt[2]))
        wire_outer.Close()
        builder.AddWire(wire_outer.Wire())
        
        builder.Build()
        
        if builder.IsDone():
            return cq.Shape(builder.Shape())
        return None
        
    except Exception:
        return None


def _create_convex_hull_hub_3strut(
    vertex: Point3D,
    radial: Point3D,
    strut_directions: List[Point3D],
    hub_inset: float,
    strut_width: float,
    strut_depth: float
) -> Optional[cq.Shape]:
    """
    Create a hub as convex hull of all strut end face corners.
    
    This guarantees the hub faces are flush with strut ends since
    the hull includes the exact corner points of each strut's end face.
    """
    half_d = strut_depth / 2.0
    half_w = strut_width / 2.0
    
    # Collect all corner points from all strut end faces
    all_corners = []
    
    for direction in strut_directions:
        # Strut endpoint (where the shortened strut ends)
        endpoint = add(vertex, scale(direction, hub_inset))
        
        # Compute strut profile corners at this endpoint
        side_normal = normalize(cross(radial, direction))
        if norm(side_normal) < 1e-9:
            if abs(radial[2]) < 0.9:
                side_normal = normalize(cross(radial, (0, 0, 1)))
            else:
                side_normal = normalize(cross(radial, (1, 0, 0)))
        
        # 4 corners of strut end face
        inner_center = sub(endpoint, scale(radial, half_d))
        outer_center = add(endpoint, scale(radial, half_d))
        
        all_corners.extend([
            sub(inner_center, scale(side_normal, half_w)),
            add(inner_center, scale(side_normal, half_w)),
            add(outer_center, scale(side_normal, half_w)),
            sub(outer_center, scale(side_normal, half_w)),
        ])
    
    # Add vertex center points for better hull shape
    all_corners.append(sub(vertex, scale(radial, half_d)))
    all_corners.append(add(vertex, scale(radial, half_d)))
    
    # Build convex hull
    try:
        from scipy.spatial import ConvexHull
        import numpy as np
        
        points = np.array(all_corners)
        hull = ConvexHull(points)
        
        # Build solid from hull faces
        sewing = BRepBuilderAPI_Sewing(1e-6)
        
        for simplex in hull.simplices:
            p0, p1, p2 = points[simplex[0]], points[simplex[1]], points[simplex[2]]
            
            poly = BRepBuilderAPI_MakePolygon()
            poly.Add(gp_Pnt(float(p0[0]), float(p0[1]), float(p0[2])))
            poly.Add(gp_Pnt(float(p1[0]), float(p1[1]), float(p1[2])))
            poly.Add(gp_Pnt(float(p2[0]), float(p2[1]), float(p2[2])))
            poly.Close()
            
            if poly.IsDone():
                face = BRepBuilderAPI_MakeFace(poly.Wire())
                if face.IsDone():
                    sewing.Add(face.Face())
        
        sewing.Perform()
        shell = sewing.SewedShape()
        
        solid_maker = BRepBuilderAPI_MakeSolid()
        solid_maker.Add(TopoDS.Shell_s(shell))
        
        if solid_maker.IsDone():
            return cq.Shape(solid_maker.Solid())
        return cq.Shape(shell)
        
    except Exception:
        return None


def _create_bar_hub(
    vertex: Point3D,
    radial: Point3D,
    strut_directions: List[Point3D],
    hub_inset: float,
    strut_width: float,
    strut_depth: float
) -> Optional[cq.Shape]:
    """
    Create a bar-shaped hub for a 2-strut vertex (typically at dome edge).
    
    This is a simple rectangular bar connecting the two strut ends.
    """
    if len(strut_directions) != 2:
        return None
    
    half_d = strut_depth / 2.0
    half_w = strut_width / 2.0
    
    d1, d2 = strut_directions[0], strut_directions[1]
    
    # Bar direction is the average of the two strut directions
    bar_dir = normalize(add(d1, d2))
    
    # If struts are opposite, use perpendicular
    if norm(bar_dir) < 1e-9:
        bar_dir = normalize(cross(radial, d1))
    
    # Side direction perpendicular to bar and radial
    side_dir = normalize(cross(radial, bar_dir))
    
    # Bar length should span from one strut end to the other
    bar_length = hub_inset * 2.0
    
    # Create 8 corners of the bar
    corners = []
    for r_sign in [-1, 1]:  # inner/outer
        for b_sign in [-1, 1]:  # along bar
            for s_sign in [-1, 1]:  # side
                corner = add(vertex,
                    add(scale(radial, r_sign * half_d),
                    add(scale(bar_dir, b_sign * bar_length * 0.4),
                        scale(side_dir, s_sign * half_w))))
                corners.append(corner)
    
    # Create box using loft between two end profiles
    try:
        # Define the two end rectangles
        end1 = [corners[0], corners[1], corners[3], corners[2]]  # inner end
        end2 = [corners[4], corners[5], corners[7], corners[6]]  # outer end
        
        builder = BRepOffsetAPI_ThruSections(True, True)
        
        wire1 = BRepBuilderAPI_MakePolygon()
        for pt in end1:
            wire1.Add(gp_Pnt(pt[0], pt[1], pt[2]))
        wire1.Close()
        builder.AddWire(wire1.Wire())
        
        wire2 = BRepBuilderAPI_MakePolygon()
        for pt in end2:
            wire2.Add(gp_Pnt(pt[0], pt[1], pt[2]))
        wire2.Close()
        builder.AddWire(wire2.Wire())
        
        builder.Build()
        
        if builder.IsDone():
            return cq.Shape(builder.Shape())
        return None
        
    except Exception:
        return None


def _create_convex_hull_hub(
    hub_info: Dict,
    strut_width: float,
    strut_depth: float,
    dome_center: Point3D
) -> Optional[cq.Shape]:
    """
    Create a convex hull hub for vertices with more than 3 struts.
    Fallback for complex cases.
    """
    vertex = hub_info['vertex']
    strut_directions = hub_info['strut_directions']
    strut_endpoints = hub_info['strut_endpoints']
    hub_inset = hub_info['hub_inset']
    
    radial = normalize(sub(vertex, dome_center))
    half_w = strut_width / 2.0
    half_d = strut_depth / 2.0
    
    all_points = []
    
    for direction, endpoint in zip(strut_directions, strut_endpoints):
        side_normal = normalize(cross(radial, direction))
        if norm(side_normal) < 1e-9:
            if abs(radial[2]) < 0.9:
                side_normal = normalize(cross(radial, (0, 0, 1)))
            else:
                side_normal = normalize(cross(radial, (1, 0, 0)))
        
        inner_center = sub(endpoint, scale(radial, half_d))
        outer_center = add(endpoint, scale(radial, half_d))
        
        all_points.extend([
            add(inner_center, scale(side_normal, half_w)),
            sub(inner_center, scale(side_normal, half_w)),
            add(outer_center, scale(side_normal, half_w)),
            sub(outer_center, scale(side_normal, half_w)),
        ])
    
    all_points.extend([
        sub(vertex, scale(radial, half_d)),
        add(vertex, scale(radial, half_d)),
    ])
    
    try:
        from scipy.spatial import ConvexHull
        import numpy as np

        points_array = np.array(all_points)
        hull = ConvexHull(points_array)
        return _build_hull_solid(points_array, hull.simplices)
    except Exception:
        return None


def _build_hull_solid(points: 'np.ndarray', simplices: 'np.ndarray') -> Optional[cq.Shape]:
    """Build a solid from convex hull triangular faces."""
    try:
        sewing = BRepBuilderAPI_Sewing(1e-6)
        
        for simplex in simplices:
            p0 = points[simplex[0]]
            p1 = points[simplex[1]]
            p2 = points[simplex[2]]
            
            poly = BRepBuilderAPI_MakePolygon()
            poly.Add(gp_Pnt(float(p0[0]), float(p0[1]), float(p0[2])))
            poly.Add(gp_Pnt(float(p1[0]), float(p1[1]), float(p1[2])))
            poly.Add(gp_Pnt(float(p2[0]), float(p2[1]), float(p2[2])))
            poly.Close()
            
            if poly.IsDone():
                face = BRepBuilderAPI_MakeFace(poly.Wire())
                if face.IsDone():
                    sewing.Add(face.Face())
        
        sewing.Perform()
        sewn_shape = sewing.SewedShape()
        
        solid_maker = BRepBuilderAPI_MakeSolid()
        solid_maker.Add(TopoDS.Shell_s(sewn_shape))
        
        if solid_maker.IsDone():
            return cq.Shape(solid_maker.Solid())
        
        return cq.Shape(sewn_shape)
        
    except Exception:
        return None


def _create_cylindrical_core_hub(
    vertex: Point3D,
    radial: Point3D,
    strut_directions: List[Point3D],
    hub_inset: float,
    strut_width: float,
    strut_depth: float,
    core_radius: Optional[float] = None
) -> Optional[cq.Shape]:
    """
    Create a cylindrical core hub for miter-cut struts.
    
    The cylinder is centered on the vertex with axis along the radial direction.
    Struts with compound miter cuts will meet this cylinder tangentially.
    
    Args:
        vertex: Hub center position
        radial: Radial direction (cylinder axis)
        strut_directions: Directions of incoming struts
        hub_inset: How far struts are shortened
        strut_width: Width of struts
        strut_depth: Depth of struts  
        core_radius: Cylinder radius (auto-computed if None)
    
    Returns:
        CadQuery Shape of the cylindrical hub
    """
    half_d = strut_depth / 2.0
    
    # Auto-compute core radius if not provided
    # The cylinder should be large enough to meet the miter-cut strut ends
    if core_radius is None:
        # For miter-cut struts meeting at a vertex, the strut end faces
        # are at approximately hub_inset distance from the vertex.
        # The cylinder radius should be close to hub_inset so the
        # strut end faces meet the cylinder surface.
        #
        # We compute the minimum distance from vertex to strut end faces
        # accounting for the strut coming in at an angle.
        min_perpendicular_dist = float('inf')
        
        for direction in strut_directions:
            # The strut end center is at hub_inset along the strut direction
            # The perpendicular distance from vertex to the strut end face
            # (which lies in a radial plane) depends on the strut angle
            
            # cos(theta) = direction . radial
            cos_theta = abs(direction[0]*radial[0] + direction[1]*radial[1] + direction[2]*radial[2])
            sin_theta = math.sqrt(max(0, 1 - cos_theta * cos_theta))
            
            # Perpendicular distance = hub_inset * sin(theta)
            if sin_theta > 0.1:
                perp_dist = hub_inset * sin_theta
                min_perpendicular_dist = min(min_perpendicular_dist, perp_dist)
        
        # Use 95% of minimum perpendicular distance for slight overlap
        if min_perpendicular_dist < float('inf'):
            core_radius = min_perpendicular_dist * 0.95
        else:
            core_radius = hub_inset * 0.8
        
        # Ensure minimum radius for structural integrity
        core_radius = max(core_radius, strut_width * 0.4)
    
    # Cylinder height spans from inner to outer face of strut depth
    cylinder_height = strut_depth
    
    try:
        # Create cylinder centered on vertex, axis along radial
        # CadQuery makes a cylinder along Z, so we need to transform it
        
        # Create cylinder
        cylinder = cq.Workplane("XY").cylinder(cylinder_height, core_radius)
        
        # We need to position and orient this cylinder
        # The cylinder should be centered on the vertex with axis along radial
        
        # Create transformation: translate to vertex and rotate to align with radial
        
        # Compute rotation to align Z with radial
        z_axis = (0, 0, 1)
        
        # Rotation axis = Z x radial
        rot_axis = cross(z_axis, radial)
        rot_axis_norm = norm(rot_axis)
        
        if rot_axis_norm > 1e-9:
            rot_axis = scale(rot_axis, 1.0 / rot_axis_norm)
            # Rotation angle = acos(Z . radial)
            dot_val = z_axis[0]*radial[0] + z_axis[1]*radial[1] + z_axis[2]*radial[2]
            rot_angle = math.acos(max(-1, min(1, dot_val)))
            rot_angle_deg = math.degrees(rot_angle)
            
            # Apply rotation and translation
            cylinder = cylinder.rotate((0, 0, 0), rot_axis, rot_angle_deg)
        
        # Translate to vertex position
        cylinder = cylinder.translate(vertex)
        
        return cylinder.val()
        
    except Exception:
        return None


def create_hub_by_style(
    hub_info: Dict,
    strut_width: float,
    strut_depth: float,
    dome_center: Point3D,
    hub_style: str = "convex_hull"
) -> Optional[cq.Shape]:
    """
    Create a hub using the specified style.
    
    Args:
        hub_info: Hub geometry from compute_hub_geometry()
        strut_width: Width of struts
        strut_depth: Depth of struts
        dome_center: Center of dome
        hub_style: One of "convex_hull", "tapered_prism", "cylindrical_core"
    
    Returns:
        CadQuery Shape for the hub, or None if invalid
    """
    vertex = hub_info['vertex']
    strut_directions = hub_info['strut_directions']
    hub_inset = hub_info['hub_inset']
    num_struts = hub_info['num_struts']
    
    if num_struts < 2:
        return None
    
    # Radial direction at the vertex
    radial = normalize(sub(vertex, dome_center))
    
    if hub_style == "tapered_prism":
        if num_struts == 3:
            result = _create_tapered_prism_hub(
                vertex, radial, strut_directions, hub_inset, strut_width, strut_depth
            )
            if result is not None:
                return result
        # Fall back to convex hull for non-3-strut or if tapered prism fails
        return create_hub_solid(hub_info, strut_width, strut_depth, dome_center)
    
    elif hub_style == "cylindrical_core":
        if num_struts >= 2:
            result = _create_cylindrical_core_hub(
                vertex, radial, strut_directions, hub_inset, strut_width, strut_depth
            )
            if result is not None:
                return result
        return create_hub_solid(hub_info, strut_width, strut_depth, dome_center)
    
    else:  # "convex_hull" or default
        return create_hub_solid(hub_info, strut_width, strut_depth, dome_center)

