"""
Selective Visualization Helpers for Geodesic Domes.

This module provides functions for generating dome parts separated by category
(complete vs. cut/partial) for visualization and manufacturing analysis.
"""

import cadquery as cq
from typing import List, Dict, Tuple, Optional

from .geometry import (
    Point3D,
    Edge,
    normalize,
    sub,
    add,
    scale,
    dot,
    build_vertex_to_edges_map,
    compute_tangent_basis,
)

from .cuboid_struts import create_shortened_cuboid_strut
from .hubs import compute_hub_geometry, create_hub_by_style, calculate_hub_corner


def generate_window_plates_by_face(
    vertices: List[Point3D],
    faces: List[List[int]],
    strut_depth: float,
    strut_width: float,
    plate_depth: float,
    dome_center: Point3D,
    margin: float,
    hub_inset: float,
) -> Dict[int, cq.Solid]:
    """
    Generate window plates indexed by face index.
    
    This function creates window plates that fit inside each window opening,
    calculated from the actual hub corner positions.
    
    Args:
        vertices: All vertex coordinates
        faces: List of faces (lists of vertex indices)
        strut_depth: Depth of struts (radial thickness)
        strut_width: Width of struts
        plate_depth: Thickness of the window plate
        dome_center: Center of the dome
        margin: Gap between plate and frame
        hub_inset: Distance struts are shortened from vertices
        
    Returns:
        Dict mapping face_idx to window plate solid (before boundary cutting)
    """
    plates_by_face = {}
    
    for face_idx, face in enumerate(faces):
        if len(face) < 3:
            continue
            
        # Calculate centroid and direction
        face_verts_coords = [vertices[i] for i in face]
        cx = sum(v[0] for v in face_verts_coords) / len(face_verts_coords)
        cy = sum(v[1] for v in face_verts_coords) / len(face_verts_coords)
        cz = sum(v[2] for v in face_verts_coords) / len(face_verts_coords)
        centroid = (cx, cy, cz)
        
        # Window direction (from origin to centroid)
        window_dir = normalize(sub(centroid, dome_center))
        
        # Calculate the actual hub corners for this window
        hub_corners = []
        n = len(face)
        
        for i in range(n):
            curr_idx = face[i]
            prev_idx = face[(i - 1 + n) % n]
            next_idx = face[(i + 1) % n]
            
            vertex = vertices[curr_idx]
            
            d_prev = normalize(sub(vertices[prev_idx], vertex))
            d_next = normalize(sub(vertices[next_idx], vertex))
            
            radial = normalize(sub(vertex, dome_center))
            
            corner = calculate_hub_corner(
                vertex, radial, d_prev, d_next, hub_inset, strut_depth
            )
            hub_corners.append(corner)
            
        # Determine window plane distance from origin
        u_values = [dot(sub(p, dome_center), window_dir) for p in hub_corners]
        u_min = min(u_values)
        
        plane_center = add(dome_center, scale(window_dir, u_min))
        
        # Tangent basis for the window plane
        u_vec, v_vec = compute_tangent_basis(window_dir)
        
        # Project corners to 2D on this plane
        points_2d = []
        for p in hub_corners:
            rel = sub(p, plane_center)
            x2d = dot(rel, u_vec)
            y2d = dot(rel, v_vec)
            points_2d.append((x2d, y2d))
            
        # Create plate
        try:
            wp = cq.Workplane(cq.Plane(origin=plane_center, xDir=u_vec, normal=window_dir))
            poly = wp.polyline(points_2d).close()
            offset_val = -margin
            plate = poly.offset2D(offset_val).extrude(plate_depth)
            plates_by_face[face_idx] = plate.val()
        except Exception:
            pass
            
    return plates_by_face


def separate_dome_parts(
    vertices: List[Point3D],
    edges: List[Edge],
    faces: List[List[int]],
    radius_cm: float,
    portion: float,
    strut_width: float,
    strut_depth: float,
    hub_inset: float,
    window_plate_depth: float,
    window_margin: float,
    hub_style: str = "tapered_prism",
    tolerance: float = 0.1,
) -> Dict:
    """
    Generate dome parts separated into complete and partial/cut categories.
    
    This function is useful for visualization and manufacturing planning,
    allowing you to see which parts are fully above the base and which
    get trimmed at the boundary.
    
    Args:
        vertices: All vertex coordinates
        edges: List of edges (vertex index pairs)
        faces: List of faces (lists of vertex indices)
        radius_cm: Dome radius in cm
        portion: Dome portion (0.5 = hemisphere)
        strut_width: Width of struts
        strut_depth: Depth of struts
        hub_inset: Distance struts are shortened from vertices
        window_plate_depth: Thickness of window plates
        window_margin: Gap between plate and frame
        hub_style: Hub style ("tapered_prism", "convex_hull", "cylindrical_core")
        tolerance: Floating point tolerance for boundary detection
    
    Returns:
        Dict with:
            - complete_struts: Compound of struts fully above base
            - partial_struts: Compound of struts cut at base
            - complete_hubs: Compound of interior hubs
            - boundary_hubs: Compound of boundary hubs (still complete pieces)
            - complete_panes: Compound of window panes fully above base
            - partial_panes: Compound of window panes cut at base
            - counts: Dict with counts for each category
    """
    dome_center = (0.0, 0.0, 0.0)
    cutoff_y = -radius_cm * (2 * portion - 1)
    cut_plane_y = cutoff_y - (strut_width / 2.0)
    
    # Create cutter for boundary trimming
    huge_dim = radius_cm * 10
    lx, ly, lz = 2 * huge_dim, cut_plane_y - (-huge_dim), 2 * huge_dim
    cutter_solid = cq.Solid.makeBox(lx, ly, lz).translate(
        cq.Vector(-huge_dim, -huge_dim, -huge_dim)
    )
    
    # Build vertex to edges map
    vertex_to_edges = build_vertex_to_edges_map(edges)
    
    # Separate struts
    complete_strut_shapes = []
    partial_strut_shapes = []
    
    for v1_idx, v2_idx in edges:
        v1, v2 = vertices[v1_idx], vertices[v2_idx]
        
        # Check if strut is complete or partial
        v1_above = v1[1] >= cutoff_y - tolerance
        v2_above = v2[1] >= cutoff_y - tolerance
        
        strut = create_shortened_cuboid_strut(
            v1, v2, strut_width, strut_depth, hub_inset, hub_inset, dome_center
        )
        
        if v1_above and v2_above:
            complete_strut_shapes.append(strut)
        else:
            # Cut the strut at the plane
            try:
                cut_strut = strut.cut(cutter_solid)
                if cut_strut.isValid() and cut_strut.Volume() > 1e-6:
                    partial_strut_shapes.append(cut_strut)
            except Exception:
                pass
    
    # Separate hubs
    complete_hub_shapes = []
    boundary_hub_shapes = []
    
    used_vertices = set()
    for v1_idx, v2_idx in edges:
        used_vertices.add(v1_idx)
        used_vertices.add(v2_idx)
    
    for v_idx in used_vertices:
        vertex = vertices[v_idx]
        if vertex[1] < cutoff_y - tolerance:
            continue
            
        # Check if this is a boundary hub
        edge_indices = vertex_to_edges.get(v_idx, [])
        neighbors_below = 0
        for edge_idx in edge_indices:
            v1_idx_e, v2_idx_e = edges[edge_idx]
            other_idx = v2_idx_e if v1_idx_e == v_idx else v1_idx_e
            if vertices[other_idx][1] < cutoff_y - tolerance:
                neighbors_below += 1
        
        hub_info = compute_hub_geometry(
            vertex, v_idx, vertices, edges, vertex_to_edges,
            strut_width, strut_depth, hub_inset, dome_center
        )
        
        hub = create_hub_by_style(hub_info, strut_width, strut_depth, dome_center, hub_style)
        if hub is not None:
            if neighbors_below > 0:
                boundary_hub_shapes.append(hub)
            else:
                complete_hub_shapes.append(hub)
    
    # Separate window panes
    plates_by_face = generate_window_plates_by_face(
        vertices, faces, strut_depth, strut_width, window_plate_depth,
        dome_center, window_margin, hub_inset
    )
    
    complete_pane_shapes = []
    partial_pane_shapes = []
    
    for face_idx, plate in plates_by_face.items():
        face = faces[face_idx]
        face_verts = [vertices[i] for i in face]
        
        # Check if all vertices are above cutoff
        all_above = all(v[1] >= cutoff_y - tolerance for v in face_verts)
        
        if all_above:
            complete_pane_shapes.append(plate)
        else:
            # Cut the pane at the plane
            try:
                cut_pane = plate.cut(cutter_solid)
                if cut_pane.isValid() and cut_pane.Volume() > 1e-6:
                    partial_pane_shapes.append(cut_pane)
            except Exception:
                pass
    
    return {
        'complete_struts': cq.Compound.makeCompound(complete_strut_shapes) if complete_strut_shapes else None,
        'partial_struts': cq.Compound.makeCompound(partial_strut_shapes) if partial_strut_shapes else None,
        'complete_hubs': cq.Compound.makeCompound(complete_hub_shapes) if complete_hub_shapes else None,
        'boundary_hubs': cq.Compound.makeCompound(boundary_hub_shapes) if boundary_hub_shapes else None,
        'complete_panes': cq.Compound.makeCompound(complete_pane_shapes) if complete_pane_shapes else None,
        'partial_panes': cq.Compound.makeCompound(partial_pane_shapes) if partial_pane_shapes else None,
        'counts': {
            'complete_struts': len(complete_strut_shapes),
            'partial_struts': len(partial_strut_shapes),
            'complete_hubs': len(complete_hub_shapes),
            'boundary_hubs': len(boundary_hub_shapes),
            'complete_panes': len(complete_pane_shapes),
            'partial_panes': len(partial_pane_shapes),
        }
    }


# Color palette for dome visualization
VISUALIZATION_COLORS = {
    'complete_struts': (101, 67, 33),       # Dark brown
    'partial_struts': (255, 140, 0),        # Dark orange
    'complete_hubs': (181, 101, 29),        # Light brown
    'boundary_hubs': (181, 101, 29),        # Light brown
    'complete_panes': (144, 238, 144),      # Light green
    'partial_panes': (255, 140, 0),         # Dark orange
}

# Alpha (transparency) values for dome visualization (0.0 = transparent, 1.0 = opaque)
VISUALIZATION_ALPHAS = {
    'complete_struts': 1.0,
    'partial_struts': 1.0,
    'complete_hubs': 1.0,
    'boundary_hubs': 1.0,
    'complete_panes': 0.5,                  # Semi-transparent
    'partial_panes': 1.0,
}

VISUALIZATION_NAMES = {
    'complete_struts': 'Complete Struts',
    'complete_hubs': 'Complete Hubs',
    'partial_struts': 'Partial Struts',
    'boundary_hubs': 'Boundary Hubs',
    'complete_panes': 'Complete Panes',
    'partial_panes': 'Partial Panes',
}

# Stained glass visualization: dark brown structure with colorful panes
STAINED_GLASS_COLORS = {
    'structure': (101, 67, 33),              # Dark brown for all structural elements
    'pane_palette': [
        (255, 215, 0),                        # Yellow (gold)
        (220, 20, 60),                        # Red (crimson)
        (30, 144, 255),                       # Blue (dodger blue)
        (50, 205, 50),                        # Green (lime green)
        (148, 0, 211),                        # Purple (dark violet)
    ],
}


def generate_stained_glass_colors(num_panes: int, seed: int = 42) -> list:
    """
    Generate a list of random colors for panes from the stained glass palette.
    
    Args:
        num_panes: Number of panes to color
        seed: Random seed for reproducibility
        
    Returns:
        List of RGB tuples
    """
    import random
    random.seed(seed)
    palette = STAINED_GLASS_COLORS['pane_palette']
    return [random.choice(palette) for _ in range(num_panes)]

