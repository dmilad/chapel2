"""
Chapel2 - CadQuery Geodesic Dome Generator

A toolkit for generating geodesic dome structures with wedged struts using CadQuery.
"""

from .geometry import (
    Point3D,
    Edge,
    Face,
    generate_geodesic_dome,
    generate_honeycomb_dome,
    analyze_polygon_sizes,
    generate_assembly_guide,
    format_assembly_guide_text,
)

from .wedged_struts import (
    create_strut_from_edge,
    create_wedged_strut,
)

from .cuboid_struts import (
    create_cuboid_strut,
    create_cuboid_strut_from_edge,
)

from .hubs import (
    compute_hub_geometry,
    create_hub_solid,
    create_hub_by_style,
)

from .miter_struts import (
    create_miter_cut_strut,
    compute_miter_cut_angles,
    analyze_miter_cuts_for_dome,
)

from .dome_generator import (
    generate_dome,
    generate_dome_with_hubs,
    generate_dome_struts_individually,
    generate_2v_test_dome,
    generate_chapel_dome,
    export_dome,
)

__version__ = "0.1.0"
__all__ = [
    # Types
    "Point3D",
    "Edge", 
    "Face",
    # Geometry functions
    "generate_geodesic_dome",
    "generate_honeycomb_dome",
    "analyze_polygon_sizes",
    "generate_assembly_guide",
    "format_assembly_guide_text",
    # Strut functions
    "create_strut_from_edge",
    "create_wedged_strut",
    "create_cuboid_strut",
    "create_cuboid_strut_from_edge",
    # Hub functions
    "compute_hub_geometry",
    "create_hub_solid",
    "create_hub_by_style",
    # Miter strut functions
    "create_miter_cut_strut",
    "compute_miter_cut_angles",
    "analyze_miter_cuts_for_dome",
    # Dome generator functions
    "generate_dome",
    "generate_dome_with_hubs",
    "generate_dome_struts_individually",
    "generate_2v_test_dome",
    "generate_chapel_dome",
    "export_dome",
]
