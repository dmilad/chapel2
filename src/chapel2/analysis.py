"""
Manufacturability Analysis for Geodesic Domes.

Provides functions to analyze unique parts, counts, and cut vs complete parts.
"""

import math
from collections import Counter
from typing import List, Dict, Tuple, Any, Set

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
    build_vertex_to_edges_map,
)


def classify_struts(
    vertices: List[Point3D],
    edges: List[Edge],
    cutoff_y: float,
    tolerance: float = 0.1,
    length_tol_cm: float = 0.1
) -> Dict[str, Any]:
    """
    Classify struts as complete or cut based on cutoff plane.
    
    Complete struts have both endpoints above the cutoff.
    Cut struts have one endpoint below the cutoff (they get trimmed).
    
    Args:
        vertices: All vertex coordinates
        edges: List of edges (vertex index pairs)
        cutoff_y: Y-coordinate cutoff (parts below this get cut)
        tolerance: Floating point tolerance for boundary detection
        length_tol_cm: Tolerance for grouping lengths
    
    Returns:
        Dict with 'complete' and 'cut' strut analysis
    """
    complete_struts = []
    cut_struts = []
    
    for edge_idx, (v1_idx, v2_idx) in enumerate(edges):
        v1 = vertices[v1_idx]
        v2 = vertices[v2_idx]
        
        # Calculate strut length
        length = norm(sub(v2, v1))
        length_rounded = round(length / length_tol_cm) * length_tol_cm
        
        # Classify based on whether both endpoints are above cutoff
        v1_above = v1[1] >= cutoff_y - tolerance
        v2_above = v2[1] >= cutoff_y - tolerance
        
        strut_info = {
            'edge_idx': edge_idx,
            'v1_idx': v1_idx,
            'v2_idx': v2_idx,
            'length_cm': length,
            'length_rounded': length_rounded,
        }
        
        if v1_above and v2_above:
            complete_struts.append(strut_info)
        else:
            # At least one endpoint is below cutoff - this strut gets cut
            cut_struts.append(strut_info)
    
    # Group by length
    complete_by_length = Counter(s['length_rounded'] for s in complete_struts)
    cut_by_length = Counter(s['length_rounded'] for s in cut_struts)
    
    return {
        'complete': {
            'count': len(complete_struts),
            'struts': complete_struts,
            'by_length': dict(complete_by_length),
            'unique_types': len(complete_by_length),
        },
        'cut': {
            'count': len(cut_struts),
            'struts': cut_struts,
            'by_length': dict(cut_by_length),
            'unique_types': len(cut_by_length),
        },
        'total': len(edges),
    }


def classify_hubs(
    vertices: List[Point3D],
    edges: List[Edge],
    vertex_to_edges: Dict[int, List[int]],
    cutoff_y: float,
    tolerance: float = 0.1,
    angle_tol_deg: float = 1.0
) -> Dict[str, Any]:
    """
    Classify hubs as interior or boundary based on position.
    
    Interior hubs: all connected struts are complete (fully above cutoff).
    Boundary hubs: at the base ring, some connected struts are trimmed.
    
    Note: Boundary hubs are still complete hub pieces - they're not cut.
    They're just located at the dome's edge where attached struts get trimmed.
    
    Hubs are characterized by number of struts and angles between them.
    
    Args:
        vertices: All vertex coordinates
        edges: List of edges (vertex index pairs)
        vertex_to_edges: Map from vertex index to edge indices
        cutoff_y: Y-coordinate cutoff
        tolerance: Floating point tolerance
        angle_tol_deg: Tolerance for grouping angles
    
    Returns:
        Dict with 'interior' and 'boundary' hub analysis
    """
    interior_hubs = []
    boundary_hubs = []
    
    for v_idx, edge_indices in vertex_to_edges.items():
        vertex = vertices[v_idx]
        num_struts = len(edge_indices)
        
        if num_struts < 2:
            continue
        
        # Get strut directions and check if neighbors are above/below cutoff
        directions = []
        neighbors_below_cutoff = 0
        
        for edge_idx in edge_indices:
            v1_idx, v2_idx = edges[edge_idx]
            other_idx = v2_idx if v1_idx == v_idx else v1_idx
            other = vertices[other_idx]
            direction = normalize(sub(other, vertex))
            directions.append(direction)
            
            if other[1] < cutoff_y - tolerance:
                neighbors_below_cutoff += 1
        
        # Compute angles between all pairs of struts
        angles = []
        for i in range(len(directions)):
            for j in range(i + 1, len(directions)):
                d1, d2 = directions[i], directions[j]
                cos_angle = dot(d1, d2)
                angle_deg = math.degrees(math.acos(max(-1, min(1, cos_angle))))
                angle_rounded = round(angle_deg / angle_tol_deg) * angle_tol_deg
                angles.append(angle_rounded)
        
        angles_sorted = tuple(sorted(angles))
        signature = (num_struts, angles_sorted)
        
        hub_info = {
            'vertex_idx': v_idx,
            'vertex': vertex,
            'num_struts': num_struts,
            'angles': angles_sorted,
            'signature': signature,
            'neighbors_below': neighbors_below_cutoff,
        }
        
        # Classify: boundary hubs have some struts that extend below cutoff
        vertex_above = vertex[1] >= cutoff_y - tolerance
        is_boundary = neighbors_below_cutoff > 0
        
        if vertex_above and not is_boundary:
            interior_hubs.append(hub_info)
        elif vertex_above:
            # Vertex is above but some neighbors are below - this is a boundary hub
            boundary_hubs.append(hub_info)
    
    # Group by signature
    interior_by_sig = Counter(h['signature'] for h in interior_hubs)
    boundary_by_sig = Counter(h['signature'] for h in boundary_hubs)
    
    return {
        'interior': {
            'count': len(interior_hubs),
            'hubs': interior_hubs,
            'by_signature': dict(interior_by_sig),
            'unique_types': len(interior_by_sig),
        },
        'boundary': {
            'count': len(boundary_hubs),
            'hubs': boundary_hubs,
            'by_signature': dict(boundary_by_sig),
            'unique_types': len(boundary_by_sig),
        },
        'total': len(interior_hubs) + len(boundary_hubs),
    }


def classify_windows(
    vertices: List[Point3D],
    faces: List[List[int]],
    cutoff_y: float,
    tolerance: float = 0.1,
    size_tol_cm: float = 0.5
) -> Dict[str, Any]:
    """
    Classify windows (faces/openings) as complete or cut based on cutoff plane.
    
    Complete windows have all vertices above the cutoff.
    Cut windows have some vertices below the cutoff (they get trimmed).
    
    Windows are characterized by number of sides and side-to-side size.
    
    Args:
        vertices: All vertex coordinates
        faces: List of faces (lists of vertex indices)
        cutoff_y: Y-coordinate cutoff
        tolerance: Floating point tolerance
        size_tol_cm: Tolerance for grouping sizes
    
    Returns:
        Dict with 'complete' and 'cut' window analysis
    """
    complete_windows = []
    cut_windows = []
    
    for face_idx, face in enumerate(faces):
        if len(face) < 3:
            continue
        
        face_verts = [vertices[i] for i in face]
        n = len(face)
        
        # Calculate centroid
        cx = sum(v[0] for v in face_verts) / n
        cy = sum(v[1] for v in face_verts) / n
        cz = sum(v[2] for v in face_verts) / n
        
        # Calculate apothem (center to edge midpoint distance)
        apothem_sum = 0
        edge_lengths = []
        
        for i in range(n):
            v1 = face_verts[i]
            v2 = face_verts[(i + 1) % n]
            
            # Edge midpoint
            mx = (v1[0] + v2[0]) / 2
            my = (v1[1] + v2[1]) / 2
            mz = (v1[2] + v2[2]) / 2
            
            # Apothem
            apothem = math.sqrt((mx - cx)**2 + (my - cy)**2 + (mz - cz)**2)
            apothem_sum += apothem
            
            # Edge length
            edge_len = math.sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + (v2[2]-v1[2])**2)
            edge_lengths.append(edge_len)
        
        avg_apothem = apothem_sum / n
        side_to_side = 2 * avg_apothem
        side_to_side_rounded = round(side_to_side / size_tol_cm) * size_tol_cm
        avg_edge = sum(edge_lengths) / n
        avg_edge_rounded = round(avg_edge / size_tol_cm) * size_tol_cm
        
        # Signature: (num_sides, side_to_side)
        signature = (n, side_to_side_rounded)
        
        window_info = {
            'face_idx': face_idx,
            'num_sides': n,
            'side_to_side_cm': side_to_side,
            'side_to_side_rounded': side_to_side_rounded,
            'avg_edge_cm': avg_edge,
            'avg_edge_rounded': avg_edge_rounded,
            'centroid': (cx, cy, cz),
            'signature': signature,
        }
        
        # Classify: window is cut if any vertex is below cutoff
        all_above = all(v[1] >= cutoff_y - tolerance for v in face_verts)
        
        if all_above:
            complete_windows.append(window_info)
        else:
            cut_windows.append(window_info)
    
    # Group by signature
    complete_by_sig = Counter(w['signature'] for w in complete_windows)
    cut_by_sig = Counter(w['signature'] for w in cut_windows)
    
    return {
        'complete': {
            'count': len(complete_windows),
            'windows': complete_windows,
            'by_signature': dict(complete_by_sig),
            'unique_types': len(complete_by_sig),
        },
        'cut': {
            'count': len(cut_windows),
            'windows': cut_windows,
            'by_signature': dict(cut_by_sig),
            'unique_types': len(cut_by_sig),
        },
        'total': len(faces),
    }


def classify_panes(
    vertices: List[Point3D],
    faces: List[List[int]],
    strut_width: float,
    hub_inset: float,
    cutoff_y: float,
    tolerance: float = 0.1,
    size_tol_cm: float = 0.5,
    margin: float = 0.2
) -> Dict[str, Any]:
    """
    Classify window panes (plates) as complete or cut based on cutoff plane.
    
    Panes are the actual window plates that fit inside the window openings.
    They're smaller than the window openings by the margin amount.
    
    Args:
        vertices: All vertex coordinates
        faces: List of faces (lists of vertex indices)
        strut_width: Width of struts (affects pane size calculation)
        hub_inset: Hub inset distance
        cutoff_y: Y-coordinate cutoff
        tolerance: Floating point tolerance
        size_tol_cm: Tolerance for grouping sizes
        margin: Gap between pane and frame
    
    Returns:
        Dict with 'complete' and 'cut' pane analysis
    """
    complete_panes = []
    cut_panes = []
    
    for face_idx, face in enumerate(faces):
        if len(face) < 3:
            continue
        
        face_verts = [vertices[i] for i in face]
        n = len(face)
        
        # Calculate centroid
        cx = sum(v[0] for v in face_verts) / n
        cy = sum(v[1] for v in face_verts) / n
        cz = sum(v[2] for v in face_verts) / n
        
        # Calculate outer apothem (center to edge midpoint distance)
        apothem_sum = 0
        edge_lengths = []
        
        for i in range(n):
            v1 = face_verts[i]
            v2 = face_verts[(i + 1) % n]
            
            # Edge midpoint
            mx = (v1[0] + v2[0]) / 2
            my = (v1[1] + v2[1]) / 2
            mz = (v1[2] + v2[2]) / 2
            
            apothem = math.sqrt((mx - cx)**2 + (my - cy)**2 + (mz - cz)**2)
            apothem_sum += apothem
            
            edge_len = math.sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + (v2[2]-v1[2])**2)
            edge_lengths.append(edge_len)
        
        avg_outer_apothem = apothem_sum / n
        avg_edge = sum(edge_lengths) / n
        
        # Pane dimensions are reduced by strut intrusion and margin
        # The hub corner is approximately hub_inset + strut_width/2 from edge center
        # Plus margin for clearance
        inset = hub_inset + margin
        
        # Inner apothem (pane fits inside this)
        inner_apothem = avg_outer_apothem - inset
        pane_side_to_side = 2 * inner_apothem
        pane_side_to_side = max(0, pane_side_to_side)  # Prevent negative
        
        # Inner edge length
        pane_edge = avg_edge - 2 * inset
        pane_edge = max(0, pane_edge)
        
        pane_side_to_side_rounded = round(pane_side_to_side / size_tol_cm) * size_tol_cm
        pane_edge_rounded = round(pane_edge / size_tol_cm) * size_tol_cm
        
        # Signature: (num_sides, pane_side_to_side)
        signature = (n, pane_side_to_side_rounded)
        
        pane_info = {
            'face_idx': face_idx,
            'num_sides': n,
            'outer_side_to_side_cm': 2 * avg_outer_apothem,
            'pane_side_to_side_cm': pane_side_to_side,
            'pane_side_to_side_rounded': pane_side_to_side_rounded,
            'outer_edge_cm': avg_edge,
            'pane_edge_cm': pane_edge,
            'pane_edge_rounded': pane_edge_rounded,
            'centroid': (cx, cy, cz),
            'signature': signature,
        }
        
        # Classify: pane is cut if any face vertex is below cutoff
        all_above = all(v[1] >= cutoff_y - tolerance for v in face_verts)
        
        if all_above:
            complete_panes.append(pane_info)
        else:
            cut_panes.append(pane_info)
    
    # Group by signature
    complete_by_sig = Counter(p['signature'] for p in complete_panes)
    cut_by_sig = Counter(p['signature'] for p in cut_panes)
    
    return {
        'complete': {
            'count': len(complete_panes),
            'panes': complete_panes,
            'by_signature': dict(complete_by_sig),
            'unique_types': len(complete_by_sig),
        },
        'cut': {
            'count': len(cut_panes),
            'panes': cut_panes,
            'by_signature': dict(cut_by_sig),
            'unique_types': len(cut_by_sig),
        },
        'total': len(faces),
    }


def full_manufacturability_analysis(
    vertices: List[Point3D],
    edges: List[Edge],
    faces: List[List[int]],
    radius_cm: float,
    portion: float,
    strut_width: float,
    hub_inset: float,
    length_tol_cm: float = 0.1,
    angle_tol_deg: float = 1.0,
    size_tol_cm: float = 0.5,
    margin: float = 0.2
) -> Dict[str, Any]:
    """
    Perform complete manufacturability analysis on a dome.
    
    Analyzes struts, hubs, windows, and panes - classifying each as
    complete or cut, and grouping by unique dimensions.
    
    Args:
        vertices: All vertex coordinates
        edges: List of edges (active edges only)
        faces: List of faces (lists of vertex indices)
        radius_cm: Dome radius
        portion: Dome portion (0.5 = hemisphere)
        strut_width: Width of struts
        hub_inset: Hub inset distance
        length_tol_cm: Tolerance for grouping strut lengths
        angle_tol_deg: Tolerance for grouping hub angles
        size_tol_cm: Tolerance for grouping window/pane sizes
        margin: Pane margin
    
    Returns:
        Dict with complete analysis for all part types
    """
    # Calculate cutoff Y coordinate
    cutoff_y = -radius_cm * (2 * portion - 1)
    tolerance = 0.1
    
    # Build vertex to edges map
    vertex_to_edges = build_vertex_to_edges_map(edges)
    
    # Identify valid vertices (above cutoff)
    valid_vertex_indices = {i for i, v in enumerate(vertices) if v[1] >= cutoff_y - tolerance}
    
    # Filter vertex_to_edges to only include valid vertices
    vertex_to_edges_filtered = {
        v: es for v, es in vertex_to_edges.items() if v in valid_vertex_indices
    }
    
    # Perform analysis
    struts = classify_struts(vertices, edges, cutoff_y, tolerance, length_tol_cm)
    hubs = classify_hubs(vertices, edges, vertex_to_edges_filtered, cutoff_y, tolerance, angle_tol_deg)
    windows = classify_windows(vertices, faces, cutoff_y, tolerance, size_tol_cm)
    panes = classify_panes(vertices, faces, strut_width, hub_inset, cutoff_y, tolerance, size_tol_cm, margin)
    
    return {
        'struts': struts,
        'hubs': hubs,
        'windows': windows,
        'panes': panes,
        'parameters': {
            'radius_cm': radius_cm,
            'portion': portion,
            'strut_width': strut_width,
            'hub_inset': hub_inset,
            'cutoff_y': cutoff_y,
        }
    }


def _consolidate_lengths_cm(by_length: Dict[float, int], tol_cm: float = 0.3) -> Dict[float, int]:
    """
    Consolidate length entries that are within tolerance of each other.
    
    Uses 0.3 cm (~1/8") tolerance by default, suitable for manufacturing.
    
    Args:
        by_length: Dict mapping length (in cm) to count
        tol_cm: Tolerance for grouping (default 0.3 cm ≈ 1/8")
    
    Returns:
        Dict mapping representative cm value to consolidated count
    """
    # Round to tolerance grid
    consolidated = {}
    for length_cm, count in by_length.items():
        # Round to nearest tol_cm
        key = round(length_cm / tol_cm) * tol_cm
        consolidated[key] = consolidated.get(key, 0) + count
    return consolidated


def _format_length(cm: float) -> str:
    """Format a length showing both inches and cm."""
    inches = cm / 2.54
    return f"{inches:.1f}\" ({cm:.1f} cm)"


def _format_size(cm: float) -> str:
    """Format a size showing both inches and cm."""
    inches = cm / 2.54
    return f"{inches:.1f}\" / {cm:.1f} cm"


def format_analysis_report(analysis: Dict[str, Any], units: str = 'cm') -> str:
    """
    Format the manufacturability analysis as a human-readable report.
    
    Shows both inches and cm for all measurements.
    
    Args:
        analysis: Output from full_manufacturability_analysis()
        units: Ignored (kept for API compatibility) - always shows both units
    
    Returns:
        Formatted string report
    """
    lines = []
    lines.append("=" * 80)
    lines.append("MANUFACTURABILITY ANALYSIS REPORT")
    lines.append("=" * 80)
    
    params = analysis['parameters']
    radius_cm = params['radius_cm']
    lines.append(f"\nDome Parameters:")
    lines.append(f"  Radius: {radius_cm/30.48:.1f} ft ({radius_cm:.1f} cm)")
    lines.append(f"  Portion: {params['portion']}")
    lines.append(f"  Strut Width: {_format_length(params['strut_width'])}")
    lines.append(f"  Hub Inset: {_format_length(params['hub_inset'])}")
    
    # Struts
    struts = analysis['struts']
    lines.append("\n" + "-" * 80)
    lines.append("STRUTS")
    lines.append("-" * 80)
    lines.append(f"\nTotal Struts: {struts['total']}")
    
    # Consolidate strut lengths for display (in cm, then show both units)
    complete_lengths = _consolidate_lengths_cm(struts['complete']['by_length'])
    cut_lengths = _consolidate_lengths_cm(struts['cut']['by_length'])
    
    lines.append(f"  Complete: {struts['complete']['count']} ({len(complete_lengths)} unique lengths)")
    lines.append(f"  Cut (at boundary): {struts['cut']['count']} ({len(cut_lengths)} unique lengths)")
    
    lines.append("\nComplete Struts by Length:")
    for length_cm, count in sorted(complete_lengths.items()):
        lines.append(f"  {_format_length(length_cm)}: {count} pieces")
    
    if cut_lengths:
        lines.append("\nCut Struts by Original Length:")
        for length_cm, count in sorted(cut_lengths.items()):
            lines.append(f"  {_format_length(length_cm)}: {count} pieces (trimmed at base)")
    
    # Hubs - all hubs are complete pieces (none cut)
    hubs = analysis['hubs']
    lines.append("\n" + "-" * 80)
    lines.append("HUBS (all complete - none cut)")
    lines.append("-" * 80)
    lines.append(f"\nTotal Hubs: {hubs['total']}")
    
    # Combine interior and boundary hub counts
    all_hub_sigs = {}
    for sig, count in hubs['interior']['by_signature'].items():
        all_hub_sigs[sig] = all_hub_sigs.get(sig, 0) + count
    for sig, count in hubs['boundary']['by_signature'].items():
        all_hub_sigs[sig] = all_hub_sigs.get(sig, 0) + count
    
    lines.append(f"  Unique types: {len(all_hub_sigs)}")
    
    lines.append("\nHubs by Configuration:")
    for sig, count in sorted(all_hub_sigs.items()):
        num_struts = sig[0]
        angles = sig[1]
        angles_str = ", ".join(f"{a:.0f}°" for a in angles)
        lines.append(f"  {num_struts}-strut [{angles_str}]: {count} pieces")
    
    # Windows - consolidate by cm value
    windows = analysis['windows']
    lines.append("\n" + "-" * 80)
    lines.append("WINDOWS (Frame Openings)")
    lines.append("-" * 80)
    lines.append(f"\nTotal Windows: {windows['total']}")
    
    # Consolidate window signatures for display (0.5 cm tolerance for sizes)
    def consolidate_shapes(by_sig: Dict, tol_cm: float = 0.5) -> Dict:
        consolidated = {}
        for sig, count in by_sig.items():
            num_sides = sig[0]
            display_size = round(sig[1] / tol_cm) * tol_cm  # Round to tolerance grid
            key = (num_sides, display_size)
            consolidated[key] = consolidated.get(key, 0) + count
        return consolidated
    
    complete_windows = consolidate_shapes(windows['complete']['by_signature'])
    cut_windows = consolidate_shapes(windows['cut']['by_signature'])
    
    lines.append(f"  Complete: {windows['complete']['count']} ({len(complete_windows)} unique sizes)")
    lines.append(f"  Cut (at boundary): {windows['cut']['count']} ({len(cut_windows)} unique sizes)")
    
    lines.append("\nComplete Windows by Size (side-to-side):")
    for (num_sides, size_cm), count in sorted(complete_windows.items()):
        shape = "Hexagon" if num_sides == 6 else "Pentagon" if num_sides == 5 else f"{num_sides}-gon"
        lines.append(f"  {shape} {_format_size(size_cm)}: {count} pieces")
    
    if cut_windows:
        lines.append("\nCut Windows by Original Size:")
        for (num_sides, size_cm), count in sorted(cut_windows.items()):
            shape = "Hexagon" if num_sides == 6 else "Pentagon" if num_sides == 5 else f"{num_sides}-gon"
            lines.append(f"  {shape} {_format_size(size_cm)}: {count} pieces (partial)")
    
    # Panes - consolidate by cm value
    panes = analysis['panes']
    lines.append("\n" + "-" * 80)
    lines.append("PANES (Window Plates)")
    lines.append("-" * 80)
    lines.append(f"\nTotal Panes: {panes['total']}")
    
    complete_panes = consolidate_shapes(panes['complete']['by_signature'])
    cut_panes = consolidate_shapes(panes['cut']['by_signature'])
    
    lines.append(f"  Complete: {panes['complete']['count']} ({len(complete_panes)} unique sizes)")
    lines.append(f"  Cut (at boundary): {panes['cut']['count']} ({len(cut_panes)} unique sizes)")
    
    lines.append("\nComplete Panes by Size (side-to-side):")
    for (num_sides, size_cm), count in sorted(complete_panes.items()):
        shape = "Hexagon" if num_sides == 6 else "Pentagon" if num_sides == 5 else f"{num_sides}-gon"
        lines.append(f"  {shape} {_format_size(size_cm)}: {count} pieces")
    
    if cut_panes:
        lines.append("\nCut Panes by Original Size:")
        for (num_sides, size_cm), count in sorted(cut_panes.items()):
            shape = "Hexagon" if num_sides == 6 else "Pentagon" if num_sides == 5 else f"{num_sides}-gon"
            lines.append(f"  {shape} {_format_size(size_cm)}: {count} pieces (trimmed)")
    
    # Summary
    lines.append("\n" + "=" * 80)
    lines.append("SUMMARY - TOTAL UNIQUE PART TYPES")
    lines.append("=" * 80)
    lines.append(f"\nUnique strut lengths:  {struts['complete']['unique_types'] + struts['cut']['unique_types']}")
    lines.append(f"Unique hub types:      {hubs['interior']['unique_types'] + hubs['boundary']['unique_types']}")
    lines.append(f"Unique window sizes:   {windows['complete']['unique_types'] + windows['cut']['unique_types']}")
    lines.append(f"Unique pane sizes:     {panes['complete']['unique_types'] + panes['cut']['unique_types']}")
    
    lines.append("\n" + "=" * 80)
    
    return "\n".join(lines)

