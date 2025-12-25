"""
CadQuery Geodesic Dome Generator with Wedged Struts.

This module generates geodesic dome structures with wedged struts - trapezoidal 
prisms whose side faces are radial planes passing through the dome center.
"""

import cadquery as cq
from typing import List, Tuple, Dict, Optional

from .geometry import (
    Point3D,
    generate_geodesic_dome,
    generate_honeycomb_dome,
    build_vertex_to_edges_map,
)

from .wedged_struts import (
    create_strut_from_edge
)

from .cuboid_struts import (
    create_cuboid_strut_from_edge,
    create_shortened_cuboid_strut
)

from .hubs import (
    compute_hub_geometry,
    create_hub_solid,
    create_hub_by_style
)

from .miter_struts import (
    create_miter_cut_strut,
    analyze_miter_cuts_for_dome
)


# =============================================================================
# DOME GENERATION
# =============================================================================

def generate_dome_with_hubs(
    radius_cm: float,
    frequency: int,
    portion: float = 0.5,
    strut_width: float = 5.08,
    strut_depth: float = 5.08,
    dome_style: str = "honeycomb",
    hub_inset: Optional[float] = None,
    hub_style: str = "convex_hull",
    strut_style: str = "cuboid",
    max_struts: int = -1
) -> Tuple[cq.Compound, cq.Compound, Dict]:
    """
    Generate a geodesic dome with hub joints (non-overlapping struts).
    
    Struts are shortened to leave room for hub connectors at each vertex.
    
    Args:
        radius_cm: Dome radius in centimeters
        frequency: Geodesic frequency (1-6, higher = more subdivisions)
        portion: Sphere portion (0.5 = hemisphere)
        strut_width: Width of strut cross-section in cm
        strut_depth: Depth of strut (radial thickness) in cm
        dome_style: "honeycomb" (hex/pent) or "triangular"
        hub_inset: How far struts are shortened from vertices (default: 0.35x strut_width)
        hub_style: "convex_hull", "tapered_prism", or "cylindrical_core"
        strut_style: "cuboid" (standard) or "miter_cut" (for cylindrical_core hubs)
        max_struts: Maximum struts to generate (-1 = all)
    
    Returns:
        Tuple of (struts_compound, hubs_compound, info dict)
    """
    dome_center = (0.0, 0.0, 0.0)
    
    # Default hub inset based on strut width
    if hub_inset is None:
        # Hub inset = how much struts are shortened
        # Should match the hub face distance so struts meet the hub faces
        hub_inset = strut_width * 0.35
    
    # Generate dome geometry based on style
    if dome_style == "honeycomb":
        vertices, faces, edges = generate_honeycomb_dome(radius_cm, frequency, portion)
    else:
        vertices, faces, edges = generate_geodesic_dome(radius_cm, frequency, portion)
    
    # Build vertex to edges map
    vertex_to_edges = build_vertex_to_edges_map(edges)
    
    # Limit struts if requested
    edges_to_use = edges[:max_struts] if max_struts > 0 else edges
    
    # Create shortened struts
    strut_shapes = []
    
    for v1_idx, v2_idx in edges_to_use:
        start = vertices[v1_idx]
        end = vertices[v2_idx]
        
        if strut_style == "miter_cut":
            strut = create_miter_cut_strut(
                start, end, strut_width, strut_depth,
                hub_inset, hub_inset, dome_center
            )
        else:
            strut = create_shortened_cuboid_strut(
                start, end, strut_width, strut_depth,
                hub_inset, hub_inset, dome_center
            )
        strut_shapes.append(strut)
    
    # Create hubs at each vertex
    hub_shapes = []
    hub_infos = []
    
    # Find which vertices are actually used by the struts we're generating
    used_vertices = set()
    for v1_idx, v2_idx in edges_to_use:
        used_vertices.add(v1_idx)
        used_vertices.add(v2_idx)
    
    for v_idx in used_vertices:
        vertex = vertices[v_idx]
        
        # Compute hub geometry
        hub_info = compute_hub_geometry(
            vertex, v_idx, vertices, edges, vertex_to_edges,
            strut_width, strut_depth, hub_inset, dome_center
        )
        hub_infos.append(hub_info)
        
        # Create hub solid using specified style
        hub = create_hub_by_style(hub_info, strut_width, strut_depth, dome_center, hub_style)
        if hub is not None:
            hub_shapes.append(hub)
    
    # Combine into compounds
    struts_compound = cq.Compound.makeCompound(strut_shapes) if strut_shapes else cq.Compound.makeCompound([])
    hubs_compound = cq.Compound.makeCompound(hub_shapes) if hub_shapes else cq.Compound.makeCompound([])
    
    info = {
        'num_vertices': len(vertices),
        'num_edges': len(edges),
        'num_faces': len(faces),
        'num_struts_generated': len(strut_shapes),
        'num_hubs_generated': len(hub_shapes),
        'vertices': vertices,
        'edges': edges,
        'faces': faces,
        'hub_infos': hub_infos,
        'radius_cm': radius_cm,
        'frequency': frequency,
        'portion': portion,
        'strut_width': strut_width,
        'strut_depth': strut_depth,
        'hub_inset': hub_inset,
        'strut_style': strut_style,
        'hub_style': hub_style,
        'joint_style': 'hub',
    }
    
    return struts_compound, hubs_compound, info


def generate_dome(
    radius_cm: float,
    frequency: int,
    portion: float = 0.5,
    strut_width: float = 5.08,
    strut_depth: float = 5.08,
    dome_style: str = "honeycomb",
    strut_style: str = "wedged",
    joint_style: str = "overlap",
    hub_inset: Optional[float] = None,
    max_struts: int = -1
) -> Tuple[cq.Compound, Dict]:
    """
    Generate a complete geodesic dome with struts.
    
    Args:
        radius_cm: Dome radius in centimeters
        frequency: Geodesic frequency (1-6, higher = more subdivisions)
        portion: Sphere portion (0.5 = hemisphere)
        strut_width: Width of strut cross-section in cm
        strut_depth: Depth of strut (radial thickness) in cm
        dome_style: "honeycomb" (hex/pent) or "triangular"
        strut_style: "wedged" (trapezoidal, radial faces) or "cuboid" (rectangular box)
        joint_style: "overlap" (struts meet at vertices) or "hub" (struts shortened for hub joints)
        hub_inset: How far struts are shortened for hub joints (default: 1.5x strut_width)
        max_struts: Maximum struts to generate (-1 = all)
    
    Returns:
        Tuple of (compound solid, info dict with vertices/edges/faces)
        For joint_style="hub", the compound contains both struts and hubs
    """
    # If hub joint style, use the dedicated function
    if joint_style == "hub":
        struts, hubs, info = generate_dome_with_hubs(
            radius_cm=radius_cm,
            frequency=frequency,
            portion=portion,
            strut_width=strut_width,
            strut_depth=strut_depth,
            dome_style=dome_style,
            hub_inset=hub_inset,
            max_struts=max_struts
        )
        # Combine struts and hubs into one compound
        all_shapes = []
        for solid in struts.Solids():
            all_shapes.append(solid)
        for solid in hubs.Solids():
            all_shapes.append(solid)
        result = cq.Compound.makeCompound(all_shapes) if all_shapes else cq.Compound.makeCompound([])
        return result, info
    
    # Standard overlap mode
    dome_center = (0.0, 0.0, 0.0)
    
    # Generate dome geometry based on style
    if dome_style == "honeycomb":
        vertices, faces, edges = generate_honeycomb_dome(radius_cm, frequency, portion)
    else:
        vertices, faces, edges = generate_geodesic_dome(radius_cm, frequency, portion)
    
    # Limit struts if requested
    if max_struts > 0:
        edges = edges[:max_struts]
    
    # Select strut creation function based on style
    if strut_style == "cuboid":
        strut_func = create_cuboid_strut_from_edge
    else:
        strut_func = create_strut_from_edge
    
    # Create all struts
    strut_shapes = []
    
    for v1_idx, v2_idx in edges:
        start = vertices[v1_idx]
        end = vertices[v2_idx]
        
        strut = strut_func(start, end, strut_width, strut_depth, dome_center)
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
        'strut_style': strut_style,
        'joint_style': joint_style,
    }
    
    return result, info


def generate_dome_struts_individually(
    radius_cm: float,
    frequency: int,
    portion: float = 0.5,
    strut_width: float = 5.08,
    strut_depth: float = 5.08,
    dome_style: str = "honeycomb",
    strut_style: str = "wedged",
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
    
    # Select strut creation function based on style
    if strut_style == "cuboid":
        strut_func = create_cuboid_strut_from_edge
    else:
        strut_func = create_strut_from_edge
    
    # Create all struts
    strut_shapes = []
    
    for v1_idx, v2_idx in edges:
        start = vertices[v1_idx]
        end = vertices[v2_idx]
        
        strut = strut_func(start, end, strut_width, strut_depth, dome_center)
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
        'strut_style': strut_style,
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
    strut_style: str = "wedged",
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
        strut_style=strut_style,
        max_struts=max_struts
    )


def generate_chapel_dome(
    radius_ft: float = 8.0,
    frequency: int = 3,
    strut_width_in: float = 2.0,
    strut_depth_in: float = 2.0,
    strut_style: str = "wedged",
    joint_style: str = "overlap",
    hub_style: str = "convex_hull",
    hub_inset_in: Optional[float] = None
) -> Tuple[cq.Compound, Dict]:
    """
    Generate the Chapel of MOOP dome with standard parameters.
    
    Args:
        radius_ft: Dome radius in feet
        frequency: Geodesic frequency (3 recommended)
        strut_width_in: Strut width in inches
        strut_depth_in: Strut depth in inches
        strut_style: "wedged", "cuboid", or "miter_cut"
        joint_style: "overlap" (struts meet at vertices) or "hub" (struts shortened for hub joints)
        hub_style: "convex_hull", "tapered_prism", or "cylindrical_core"
        hub_inset_in: How far struts are shortened for hubs, in inches (default auto-calculated)
    
    Returns:
        Tuple of (compound, info dict)
    """
    # Convert to centimeters
    radius_cm = radius_ft * 30.48
    strut_width = strut_width_in * 2.54
    strut_depth = strut_depth_in * 2.54
    hub_inset = hub_inset_in * 2.54 if hub_inset_in is not None else None
    
    if joint_style == "hub":
        struts, hubs, info = generate_dome_with_hubs(
            radius_cm=radius_cm,
            frequency=frequency,
            portion=0.5,
            strut_width=strut_width,
            strut_depth=strut_depth,
            dome_style="honeycomb",
            hub_inset=hub_inset,
            hub_style=hub_style,
            strut_style=strut_style if strut_style in ["cuboid", "miter_cut"] else "cuboid"
        )
        # Combine struts and hubs into one compound
        all_shapes = []
        for solid in struts.Solids():
            all_shapes.append(solid)
        for solid in hubs.Solids():
            all_shapes.append(solid)
        result = cq.Compound.makeCompound(all_shapes) if all_shapes else cq.Compound.makeCompound([])
        return result, info
    
    return generate_dome(
        radius_cm=radius_cm,
        frequency=frequency,
        portion=0.5,
        strut_width=strut_width,
        strut_depth=strut_depth,
        dome_style="honeycomb",
        strut_style=strut_style,
        joint_style=joint_style,
        hub_inset=hub_inset
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
