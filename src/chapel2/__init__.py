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

from .dome_generator import (
    generate_dome,
    generate_dome_struts_individually,
    generate_2v_test_dome,
    generate_chapel_dome,
    create_strut_from_edge,
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
    # Dome generator functions
    "generate_dome",
    "generate_dome_struts_individually",
    "generate_2v_test_dome",
    "generate_chapel_dome",
    "create_strut_from_edge",
    "export_dome",
]

