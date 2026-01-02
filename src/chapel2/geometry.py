"""
Geodesic geometry calculations.
Pure Python - no Fusion 360 dependencies, so it can be tested standalone.

Ported from Chapel project for use with CadQuery.
"""

import math
from typing import List, Tuple, Dict, Set, Optional

# Type aliases for clarity
Point3D = Tuple[float, float, float]
Face = Tuple[int, int, int]
Edge = Tuple[int, int]


# =============================================================================
# VECTOR MATH HELPERS
# =============================================================================

def dot(a: Point3D, b: Point3D) -> float:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def cross(a: Point3D, b: Point3D) -> Point3D:
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])

def sub(a: Point3D, b: Point3D) -> Point3D:
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def add(a: Point3D, b: Point3D) -> Point3D:
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def scale(a: Point3D, s: float) -> Point3D:
    return (a[0]*s, a[1]*s, a[2]*s)

def norm(a: Point3D) -> float:
    return math.sqrt(dot(a, a))

def normalize(a: Point3D) -> Point3D:
    l = norm(a)
    if l < 1e-9: return (0.0, 0.0, 0.0)
    return scale(a, 1.0/l)


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


def create_icosahedron() -> Tuple[List[Point3D], List[Face]]:
    """
    Create the 12 vertices and 20 faces of a unit icosahedron.

    Returns:
        Tuple of (vertices, faces) where vertices are normalized to unit sphere
    """
    phi = (1 + math.sqrt(5)) / 2  # Golden ratio ≈ 1.618

    # Icosahedron vertices (before normalization)
    raw_vertices = [
        (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
        (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
        (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1)
    ]

    # Normalize to unit sphere
    vertices = [normalize_point(v) for v in raw_vertices]

    # 20 triangular faces (vertex indices, wound counter-clockwise)
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)
    ]

    return vertices, faces


def normalize_point(p: Point3D) -> Point3D:
    """Normalize a point to lie on the unit sphere."""
    x, y, z = p
    length = math.sqrt(x*x + y*y + z*z)
    if length == 0:
        return (0, 0, 0)
    return (x/length, y/length, z/length)


def subdivide_icosahedron(vertices: List[Point3D],
                          faces: List[Face],
                          frequency: int) -> Tuple[List[Point3D], List[Face]]:
    """
    Subdivide triangular faces to achieve desired geodesic frequency.

    Each subdivision step divides each triangle into 4 smaller triangles.

    Args:
        vertices: List of vertex coordinates
        faces: List of triangular faces (vertex index tuples)
        frequency: Target frequency (1 = no subdivision, 2 = one level, etc.)

    Returns:
        Tuple of (new_vertices, new_faces)
    """
    vertices = list(vertices)  # Make a mutable copy

    for _ in range(frequency - 1):
        new_faces = []
        edge_midpoints: Dict[Edge, int] = {}

        for face in faces:
            v0, v1, v2 = face

            # Get or create midpoints for each edge
            m01 = _get_or_create_midpoint(vertices, edge_midpoints, v0, v1)
            m12 = _get_or_create_midpoint(vertices, edge_midpoints, v1, v2)
            m20 = _get_or_create_midpoint(vertices, edge_midpoints, v2, v0)

            # Create 4 new triangles from the subdivided face
            #
            #        v0
            #       /  \
            #     m01--m20
            #     / \  / \
            #   v1--m12--v2
            #
            new_faces.extend([
                (v0, m01, m20),
                (m01, v1, m12),
                (m20, m12, v2),
                (m01, m12, m20),  # Center triangle
            ])

        faces = new_faces

    return vertices, faces


def _get_or_create_midpoint(vertices: List[Point3D],
                            cache: Dict[Edge, int],
                            v1: int,
                            v2: int) -> int:
    """
    Get existing midpoint index or create a new one.

    The midpoint is normalized to the unit sphere.
    """
    # Use sorted tuple as key so (a,b) and (b,a) map to same midpoint
    key = (min(v1, v2), max(v1, v2))

    if key not in cache:
        p1, p2 = vertices[v1], vertices[v2]
        midpoint = (
            (p1[0] + p2[0]) / 2,
            (p1[1] + p2[1]) / 2,
            (p1[2] + p2[2]) / 2,
        )
        # Normalize to unit sphere (this is what makes it geodesic)
        midpoint = normalize_point(midpoint)
        cache[key] = len(vertices)
        vertices.append(midpoint)

    return cache[key]


def scale_to_radius(vertices: List[Point3D], radius: float) -> List[Point3D]:
    """Scale unit sphere vertices to desired radius."""
    return [(x * radius, y * radius, z * radius) for x, y, z in vertices]


def filter_to_dome(vertices: List[Point3D],
                   faces: List[Face],
                   portion: float) -> Tuple[List[Point3D], List[Face]]:
    """
    Filter to keep only the upper portion of the sphere.

    Args:
        vertices: Sphere vertices
        faces: Triangular faces
        portion: What fraction of sphere to keep (0.5 = hemisphere)

    Returns:
        Tuple of (filtered_vertices, filtered_faces) with remapped indices
    """
    if not vertices:
        return [], []

    # Calculate the minimum Y value (using Y as up axis)
    # For a unit sphere: portion=0.5 means min_y=0, portion=1.0 means min_y=-1
    radius = math.sqrt(sum(c*c for c in vertices[0]))
    min_y = -radius * (2 * portion - 1)

    # Find vertices above the cutoff
    valid_vertices: Set[int] = set()
    for i, (x, y, z) in enumerate(vertices):
        if y >= min_y - 0.001:  # Small tolerance for floating point
            valid_vertices.add(i)

    # Keep faces where ALL vertices are valid
    valid_faces = [f for f in faces if all(v in valid_vertices for v in f)]

    # Remap vertex indices to remove gaps
    old_to_new: Dict[int, int] = {}
    new_vertices: List[Point3D] = []

    for face in valid_faces:
        for v in face:
            if v not in old_to_new:
                old_to_new[v] = len(new_vertices)
                new_vertices.append(vertices[v])

    new_faces = [tuple(old_to_new[v] for v in f) for f in valid_faces]

    return new_vertices, new_faces


def extract_edges(faces: List[Face]) -> List[Edge]:
    """
    Extract unique edges from faces.

    Returns edges as sorted tuples (smaller_index, larger_index).
    """
    edges: Set[Edge] = set()

    for face in faces:
        for i in range(3):
            v1, v2 = face[i], face[(i + 1) % 3]
            edge = (min(v1, v2), max(v1, v2))
            edges.add(edge)

    return list(edges)


def calculate_edge_length(vertices: List[Point3D], edge: Edge) -> float:
    """Calculate the length of an edge."""
    p1 = vertices[edge[0]]
    p2 = vertices[edge[1]]

    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    dz = p2[2] - p1[2]

    return math.sqrt(dx*dx + dy*dy + dz*dz)


def generate_cut_list(vertices: List[Point3D],
                      edges: List[Edge],
                      precision: int = 1) -> Dict[float, int]:
    """
    Generate a cut list grouping struts by length.

    Args:
        vertices: Vertex coordinates
        edges: List of edges
        precision: Decimal places to round to for grouping

    Returns:
        Dict mapping length (rounded) to count
    """
    lengths: Dict[float, int] = {}

    for edge in edges:
        length = calculate_edge_length(vertices, edge)
        length_rounded = round(length, precision)
        lengths[length_rounded] = lengths.get(length_rounded, 0) + 1

    return dict(sorted(lengths.items()))


def calculate_polygon_size(vertices: List[Point3D], face: List[int]) -> float:
    """
    Calculate the side-to-side (flat-to-flat) distance of a polygon.
    This is the distance between the midpoints of opposite edges.
    For odd-sided polygons (pentagons), it's the distance from edge midpoint to opposite vertex.

    Returns the size in the same units as the vertices.
    """
    n = len(face)
    if n < 3:
        return 0.0

    # Get polygon vertices
    poly_verts = [vertices[i] for i in face]

    # Calculate centroid
    cx = sum(v[0] for v in poly_verts) / n
    cy = sum(v[1] for v in poly_verts) / n
    cz = sum(v[2] for v in poly_verts) / n

    # Calculate average distance from centroid to edge midpoints (apothem)
    apothem_sum = 0.0
    for i in range(n):
        # Edge midpoint
        v1 = poly_verts[i]
        v2 = poly_verts[(i + 1) % n]
        mx = (v1[0] + v2[0]) / 2
        my = (v1[1] + v2[1]) / 2
        mz = (v1[2] + v2[2]) / 2

        # Distance from centroid to midpoint
        dist = math.sqrt((mx - cx)**2 + (my - cy)**2 + (mz - cz)**2)
        apothem_sum += dist

    apothem = apothem_sum / n

    # Side-to-side = 2 * apothem
    return 2 * apothem


def analyze_polygon_sizes(vertices: List[Point3D],
                          faces: List[List[int]]) -> Dict[str, Dict[str, float]]:
    """
    Analyze the sizes of hexagons and pentagons in the dome.

    Returns:
        Dict with 'hexagons' and 'pentagons' keys, each containing
        'min', 'max', 'avg' side-to-side distances
    """
    hex_sizes = []
    pent_sizes = []

    for face in faces:
        size = calculate_polygon_size(vertices, face)
        if len(face) == 6:
            hex_sizes.append(size)
        elif len(face) == 5:
            pent_sizes.append(size)

    result = {}

    if hex_sizes:
        result['hexagons'] = {
            'count': len(hex_sizes),
            'min': min(hex_sizes),
            'max': max(hex_sizes),
            'avg': sum(hex_sizes) / len(hex_sizes)
        }

    if pent_sizes:
        result['pentagons'] = {
            'count': len(pent_sizes),
            'min': min(pent_sizes),
            'max': max(pent_sizes),
            'avg': sum(pent_sizes) / len(pent_sizes)
        }

    return result


def calculate_radius_for_hexagon_size(target_hex_size: float,
                                       frequency: int,
                                       portion: float = 0.5) -> float:
    """
    Calculate the dome radius needed to achieve a target hexagon side-to-side size.

    Args:
        target_hex_size: Desired hexagon side-to-side distance
        frequency: Geodesic frequency
        portion: Dome portion

    Returns:
        Required dome radius
    """
    # Generate a unit-radius dome to measure the ratio
    verts, faces, _ = generate_honeycomb_dome(radius=1.0, frequency=frequency, portion=portion)

    sizes = analyze_polygon_sizes(verts, faces)

    if 'hexagons' not in sizes:
        return target_hex_size  # Fallback

    unit_hex_size = sizes['hexagons']['avg']

    # Scale factor to achieve target size
    scale = target_hex_size / unit_hex_size

    return scale  # This is the radius needed


def compute_dual(vertices: List[Point3D],
                 faces: List[Face]) -> Tuple[List[Point3D], List[List[int]]]:
    """
    Compute the dual of the triangular mesh (honeycomb structure).

    The dual vertices are the centroids of the original faces.
    The dual faces (polygons) correspond to the original vertices.

    Returns:
        Tuple of (dual_vertices, dual_faces)
        where dual_faces is a list of lists of vertex indices (polygons)
    """
    # 1. Compute centroids of all faces (these are the NEW vertices)
    dual_vertices = []
    for face in faces:
        v0, v1, v2 = face
        p0, p1, p2 = vertices[v0], vertices[v1], vertices[v2]
        centroid = (
            (p0[0] + p1[0] + p2[0]) / 3,
            (p0[1] + p1[1] + p2[1]) / 3,
            (p0[2] + p1[2] + p2[2]) / 3
        )
        dual_vertices.append(centroid)

    # 2. Build adjacency map: vertex_index -> list of face_indices that touch it
    # This helps us find the faces surrounding each original vertex
    vertex_to_faces: Dict[int, List[int]] = {}
    for face_idx, face in enumerate(faces):
        for v in face:
            if v not in vertex_to_faces:
                vertex_to_faces[v] = []
            vertex_to_faces[v].append(face_idx)

    # 3. Construct dual faces (polygons)
    # Each original vertex becomes a dual face (hexagon or pentagon)
    dual_faces = []

    for v_idx, touching_face_indices in vertex_to_faces.items():
        # We need to order the centroids to form a proper polygon
        # For a vertex on a sphere, the faces around it form a cycle.

        if not touching_face_indices:
            continue

        current_faces = touching_face_indices
        ordered_indices = [current_faces[0]]
        current_faces = set(current_faces[1:])

        while current_faces:
            last_face_idx = ordered_indices[-1]
            last_face = faces[last_face_idx]

            # Find a neighbor face that shares an edge connected to v_idx
            found_next = False
            for candidate_idx in list(current_faces):
                candidate_face = faces[candidate_idx]

                # Check if they share 2 vertices (common edge)
                shared = 0
                for fv in last_face:
                    if fv in candidate_face:
                        shared += 1

                # They must share 2 vertices to be neighbors
                if shared == 2:
                    ordered_indices.append(candidate_idx)
                    current_faces.remove(candidate_idx)
                    found_next = True
                    break

            if not found_next:
                # This happens for boundary vertices on a partial dome
                # where the cycle is not closed.
                # Just take what we have or break.
                # For now, we'll just append remaining (might be unordered)
                # But typically we want valid polygons.
                break

        dual_faces.append(ordered_indices)

    return dual_vertices, dual_faces


def filter_dual_to_dome(dual_vertices: List[Point3D],
                        dual_faces: List[List[int]],
                        portion: float,
                        strut_width: float = 0.0) -> Tuple[List[Point3D], List[List[int]]]:
    """
    Filter dual faces (honeycomb) to keep only the upper portion.
    
    Args:
        dual_vertices: Vertices of the dual mesh
        dual_faces: Polygon faces as lists of vertex indices
        portion: Sphere portion (0.5 = hemisphere)
        strut_width: Width of struts (currently unused, kept for API compatibility)
    
    Returns:
        Filtered (vertices, faces)
    """
    if not dual_vertices:
        return [], []

    # Approximate radius from first vertex
    radius = math.sqrt(sum(c*c for c in dual_vertices[0]))
    
    # Cutoff for the portion
    min_y = -radius * (2 * portion - 1)

    # Identify valid faces where ALL vertices are above the threshold
    valid_face_indices = []

    for i, face in enumerate(dual_faces):
        # Check if ALL vertices of this polygon are above the threshold
        all_above = all(dual_vertices[v][1] >= min_y for v in face)
        
        if all_above:
            valid_face_indices.append(i)

    # 2. Collect used vertices
    used_vertex_indices = set()
    for f_idx in valid_face_indices:
        for v in dual_faces[f_idx]:
            used_vertex_indices.add(v)

    # 3. Remap vertices
    old_to_new = {}
    new_vertices = []

    for v_idx in sorted(used_vertex_indices):
        old_to_new[v_idx] = len(new_vertices)
        new_vertices.append(dual_vertices[v_idx])

    # 4. Construct new faces
    new_faces = []
    for f_idx in valid_face_indices:
        new_face = [old_to_new[v] for v in dual_faces[f_idx]]
        new_faces.append(new_face)

    return new_vertices, new_faces


def extract_polygon_edges(faces: List[List[int]]) -> List[Edge]:
    """Extract unique edges from polygon faces."""
    edges: Set[Edge] = set()

    for face in faces:
        count = len(face)
        for i in range(count):
            v1, v2 = face[i], face[(i + 1) % count]
            edge = (min(v1, v2), max(v1, v2))
            edges.add(edge)

    return list(edges)


def generate_honeycomb_dome(radius: float,
                            frequency: int,
                            portion: float = 0.5,
                            strut_width: float = 0.0) -> Tuple[List[Point3D], List[List[int]], List[Edge]]:
    """
    Generate a honeycomb (dual) geodesic dome.

    Args:
        radius: Dome radius
        frequency: Geodesic frequency
        portion: Sphere portion (0.5 = hemisphere)
        strut_width: Strut width for boundary margin calculation (ensures struts
                     at the boundary aren't cut along their length)

    Returns:
        (vertices, faces, edges)
        faces are lists of vertex indices (polygons)
    """
    # 1. Generate base triangular mesh (Class I)
    # Note: We generate a full sphere first, then dualize, then filter
    # This ensures boundary polygons are formed correctly before cutting
    tri_vertices, tri_faces = create_icosahedron()
    if frequency > 1:
        tri_vertices, tri_faces = subdivide_icosahedron(tri_vertices, tri_faces, frequency)

    # Project to sphere
    tri_vertices = [normalize_point(v) for v in tri_vertices]

    # 2. Compute Dual
    dual_verts, dual_faces = compute_dual(tri_vertices, tri_faces)

    # 3. Scale
    dual_verts = scale_to_radius(dual_verts, radius)

    # 4. Filter with strut width margin
    final_verts, final_faces = filter_dual_to_dome(dual_verts, dual_faces, portion, strut_width)

    # 5. Edges
    final_edges = extract_polygon_edges(final_faces)

    return final_verts, final_faces, final_edges


def generate_geodesic_dome(radius: float,
                           frequency: int,
                           portion: float = 0.5) -> Tuple[List[Point3D], List[Face], List[Edge]]:
    """
    High-level function to generate a complete geodesic dome (Triangular).

    Args:
        radius: Dome radius in desired units
        frequency: Geodesic frequency (1-5 typical)
        portion: Sphere portion (0.5 = hemisphere)

    Returns:
        Tuple of (vertices, faces, edges)
    """
    # Start with icosahedron
    vertices, faces = create_icosahedron()

    # Subdivide to desired frequency
    if frequency > 1:
        vertices, faces = subdivide_icosahedron(vertices, faces, frequency)

    # Scale to desired radius
    vertices = scale_to_radius(vertices, radius)

    # Filter to dome portion
    vertices, faces = filter_to_dome(vertices, faces, portion)

    # Extract edges
    edges = extract_edges(faces)

    return vertices, faces, edges


def build_vertex_to_edges_map(edges: List[Edge]) -> Dict[int, List[int]]:
    """
    Build a mapping from each vertex index to the indices of edges connected to it.

    Args:
        edges: List of edges (vertex index pairs)

    Returns:
        Dict mapping vertex_index -> list of edge_indices
    """
    vertex_to_edges: Dict[int, List[int]] = {}

    for edge_idx, (v1, v2) in enumerate(edges):
        if v1 not in vertex_to_edges:
            vertex_to_edges[v1] = []
        if v2 not in vertex_to_edges:
            vertex_to_edges[v2] = []

        vertex_to_edges[v1].append(edge_idx)
        vertex_to_edges[v2].append(edge_idx)

    return vertex_to_edges


# =============================================================================
# WINDOW SIZE ANALYSIS FUNCTIONS
# =============================================================================

def calculate_window_dimensions(vertices: List[Point3D],
                                 face_indices: List[int],
                                 strut_width: float) -> Dict[str, float]:
    """
    Calculate outer (frame) and inner (window opening) dimensions of a polygon.

    Args:
        vertices: All vertex coordinates
        face_indices: Indices of vertices forming this face
        strut_width: Width of the struts (cm)

    Returns:
        Dict with outer/inner side-to-side and edge lengths
    """
    poly_verts = [vertices[i] for i in face_indices]
    n = len(poly_verts)

    # Calculate centroid
    cx = sum(v[0] for v in poly_verts) / n
    cy = sum(v[1] for v in poly_verts) / n
    cz = sum(v[2] for v in poly_verts) / n

    # Calculate apothem (center to edge midpoint) and edge lengths
    apothem_sum = 0
    edge_lengths = []

    for i in range(n):
        v1 = poly_verts[i]
        v2 = poly_verts[(i + 1) % n]

        # Edge midpoint
        mx = (v1[0] + v2[0]) / 2
        my = (v1[1] + v2[1]) / 2
        mz = (v1[2] + v2[2]) / 2

        # Apothem
        apothem = math.sqrt((mx-cx)**2 + (my-cy)**2 + (mz-cz)**2)
        apothem_sum += apothem

        # Edge length
        edge_len = math.sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + (v2[2]-v1[2])**2)
        edge_lengths.append(edge_len)

    avg_apothem = apothem_sum / n
    avg_edge = sum(edge_lengths) / n

    # Outer dimensions
    outer_side_to_side = 2 * avg_apothem

    # Inner window opening (strut is centered on edge)
    inner_apothem = avg_apothem - (strut_width / 2)
    inner_side_to_side = 2 * inner_apothem
    inner_edge = avg_edge - strut_width

    return {
        'sides': n,
        'outer_side_to_side': outer_side_to_side,
        'outer_edge_length': avg_edge,
        'inner_side_to_side': inner_side_to_side,
        'inner_edge_length': max(0, inner_edge),  # Prevent negative
        'face_indices': face_indices
    }


def analyze_window_sizes(vertices: List[Point3D],
                          faces: List[List[int]],
                          strut_width: float) -> Dict[str, any]:
    """
    Analyze all window opening sizes in the dome.

    Returns:
        Dictionary with hexagon and pentagon window statistics
    """
    hex_windows = []
    pent_windows = []

    for i, face in enumerate(faces):
        dims = calculate_window_dimensions(vertices, face, strut_width)
        dims['face_index'] = i

        if len(face) == 6:
            hex_windows.append(dims)
        elif len(face) == 5:
            pent_windows.append(dims)

    result = {}

    if hex_windows:
        inner_sizes = [w['inner_side_to_side'] for w in hex_windows]
        result['hexagons'] = {
            'count': len(hex_windows),
            'inner_side_to_side_min': min(inner_sizes),
            'inner_side_to_side_max': max(inner_sizes),
            'inner_side_to_side_avg': sum(inner_sizes) / len(inner_sizes),
            'outer_side_to_side_avg': sum(w['outer_side_to_side'] for w in hex_windows) / len(hex_windows),
            'windows': hex_windows
        }

    if pent_windows:
        inner_sizes = [w['inner_side_to_side'] for w in pent_windows]
        result['pentagons'] = {
            'count': len(pent_windows),
            'inner_side_to_side_min': min(inner_sizes),
            'inner_side_to_side_max': max(inner_sizes),
            'inner_side_to_side_avg': sum(inner_sizes) / len(inner_sizes),
            'outer_side_to_side_avg': sum(w['outer_side_to_side'] for w in pent_windows) / len(pent_windows),
            'windows': pent_windows
        }

    return result


# =============================================================================
# ASSEMBLY GUIDE GENERATION
# =============================================================================

def generate_assembly_guide(vertices: List[Point3D],
                            edges: List[Edge],
                            faces: List[List[int]],
                            strut_width: float) -> Dict[str, any]:
    """
    Generate a comprehensive assembly guide for the dome.

    Returns:
        Dictionary containing all assembly information
    """
    # Window analysis
    window_analysis = analyze_window_sizes(vertices, faces, strut_width)

    # Strut lengths
    strut_lengths = []
    for i, (v1, v2) in enumerate(edges):
        p1, p2 = vertices[v1], vertices[v2]
        full_length = math.sqrt(sum((a-b)**2 for a,b in zip(p1, p2)))
        # Wedged struts meet at the vertices, so cut length = full length (chord length)
        cut_length = full_length

        strut_lengths.append({
            'edge_index': i,
            'full_length': full_length,
            'cut_length': cut_length,
            'from_vertex': v1,
            'to_vertex': v2
        })

    # Group struts by cut length
    length_groups = {}
    for strut in strut_lengths:
        rounded_len = round(strut['cut_length'], 1)
        if rounded_len not in length_groups:
            length_groups[rounded_len] = []
        length_groups[rounded_len].append(strut)

    return {
        'windows': window_analysis,
        'struts': strut_lengths,
        'strut_cut_list': {k: len(v) for k, v in sorted(length_groups.items())},
        'totals': {
            'vertices': len(vertices),
            'edges': len(edges),
            'faces': len(faces)
        }
    }


def format_assembly_guide_text(guide: Dict[str, any],
                                strut_width: float,
                                strut_depth: float) -> str:
    """
    Format the assembly guide as human-readable text.
    """
    lines = []
    lines.append("=" * 70)
    lines.append("CHAPEL OF MOOP - ASSEMBLY GUIDE")
    lines.append("=" * 70)
    lines.append("")

    # Summary
    totals = guide['totals']
    lines.append("STRUCTURE SUMMARY")
    lines.append("-" * 70)
    lines.append(f"  Total Vertices:     {totals['vertices']}")
    lines.append(f"  Total Struts:       {totals['edges']}")
    lines.append(f"  Total Windows:      {totals['faces']}")
    lines.append("")

    # Strut dimensions
    lines.append("STRUT DIMENSIONS")
    lines.append("-" * 70)
    lines.append(f"  Cross-section:      {strut_width:.1f} x {strut_depth:.1f} cm ({strut_width/2.54:.1f}\" x {strut_depth/2.54:.1f}\")")
    lines.append(f"  Type:               Wedged Struts (Radial Side Faces)")
    lines.append("")

    # Window sizes
    lines.append("WINDOW OPENING SIZES (Resin Panel Dimensions)")
    lines.append("-" * 70)
    if 'hexagons' in guide['windows']:
        h = guide['windows']['hexagons']
        lines.append(f"  Hexagons ({h['count']} total):")
        lines.append(f"    Side-to-side: {h['inner_side_to_side_min']:.1f} - {h['inner_side_to_side_max']:.1f} cm")
        lines.append(f"                  ({h['inner_side_to_side_min']/2.54:.1f}\" - {h['inner_side_to_side_max']/2.54:.1f}\")")
        lines.append(f"    Average:      {h['inner_side_to_side_avg']:.1f} cm ({h['inner_side_to_side_avg']/2.54:.1f}\")")
    if 'pentagons' in guide['windows']:
        p = guide['windows']['pentagons']
        lines.append(f"  Pentagons ({p['count']} total):")
        lines.append(f"    Side-to-side: {p['inner_side_to_side_avg']:.1f} cm ({p['inner_side_to_side_avg']/2.54:.1f}\")")
    lines.append("")

    # Cut list
    lines.append("STRUT CUT LIST (Actual Cut Lengths)")
    lines.append("-" * 70)
    for length, count in guide['strut_cut_list'].items():
        lines.append(f"  {length:.1f} cm ({length/2.54:.1f}\"): {count} pieces")
    lines.append("")

    lines.append("=" * 70)

    return "\n".join(lines)
