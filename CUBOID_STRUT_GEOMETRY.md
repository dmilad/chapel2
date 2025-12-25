# Cuboid Struts and Joint Geometry in Geodesic Domes

This document explains the mathematical constraints and possibilities for using cuboid (rectangular box) struts in geodesic dome construction.

## The Problem

In a geodesic dome, multiple struts meet at each vertex. For a honeycomb (dual) geodesic dome:
- **3 struts** meet at each interior vertex
- **2 struts** meet at base/edge vertices
- The angles between struts **vary** across the dome due to geodesic geometry

### Wedged Struts (The Original Solution)

The original design uses **wedged struts** whose side faces are *radial planes* — they pass through the dome center.

```
         Dome Center
              ●
             /|\
            / | \
           /  |  \
          /   |   \     ← Radial planes converge at center
         /    |    \
        ▼     ▼     ▼
      Strut  Strut  Strut
```

**Why this works:** At any vertex, all struts' side faces converge toward the dome center. They naturally meet without gaps because they share a common focal point. No separate hub connector is needed.

### Cuboid Struts (The Challenge)

Cuboid struts have **parallel side faces** (constant rectangular cross-section).

```
    ┌────────────────┐
    │                │   ← Side faces are parallel
    │    CUBOID      │
    │                │
    └────────────────┘
```

**The problem:** When 3 cuboid struts meet at a vertex from different 3D angles, their side faces do NOT align. There's no geometric configuration where they meet flush without:
- Overlapping (physically impossible)
- Having gaps (requiring a hub/connector)

---

## Why Simple Prismatic Hubs Don't Work

### The Intuition

It seems like a triangular prism should work as a hub:

```
        Strut A
           │
      ─────┼─────
     /     │     \
    /      │      \
   ────────┼────────  Strut B
           │
      ─────┼─────
           │
        Strut C
```

But this only works if all three strut directions lie in the same plane (are **coplanar**).

### The Math

For a triangular prism hub:
- Face 1 is perpendicular to strut direction **D₁**
- Face 2 is perpendicular to strut direction **D₂**
- Face 3 is perpendicular to strut direction **D₃**

The edges of the prism are where pairs of faces intersect:
- Edge 1 = **D₁ × D₂** (cross product)
- Edge 2 = **D₂ × D₃**
- Edge 3 = **D₃ × D₁**

For this to be a **prism** (all edges parallel), we need:

```
D₁ × D₂  ∥  D₂ × D₃  ∥  D₃ × D₁
```

This condition is satisfied **only if D₁, D₂, D₃ are coplanar**.

### In a Geodesic Dome

At most vertices in a geodesic dome, the 3 strut directions are **NOT coplanar**. They spread out into 3D space because:
- The vertex sits on a curved surface
- Each strut reaches toward a different neighboring vertex
- These neighbors are at different positions on the sphere

**Conclusion:** A simple triangular prism cannot serve as a universal hub for cuboid struts in a geodesic dome.

---

## The Inner Face Triangle Solution

### Key Geometric Insight

At each vertex, there exists a **tangent plane** perpendicular to the radial direction (the local "horizontal" of the dome surface).

```
                    Radial direction
                          ↑
                          │
        ──────────────────┼──────────────────  ← Tangent plane
                          │
                       Vertex
```

If we orient each strut so its **inner face** (facing the dome center) is perpendicular to the radial direction at the vertex, then:

1. All 3 struts' inner faces are parallel to the same tangent plane
2. Their inner edges lie in the same plane
3. These edges **form a triangle**

### Visualization

Looking from inside the dome toward a vertex:

```
              Inner face of Strut 1
                    ┌─────┐
                   /       \
                  /         \
     Inner face  /           \  Inner face
     of Strut 2 /             \ of Strut 3
               └───────────────┘
               
        The 3 inner edges form a TRIANGLE
```

### The Width Constraint

The triangle is formed by the **infinite extensions** of the inner edges. Whether the actual finite-width strut edges meet depends on:

1. **Strut width (w)** — wider struts have longer edges
2. **Angles between struts (θ)** — determines edge spacing

For edges to meet perfectly, the strut width must satisfy a specific relationship with the strut angles.

#### Computing the Required Width

At a vertex with 3 strut directions **D₁, D₂, D₃** (projected onto the tangent plane as **d₁, d₂, d₃**):

The angle between struts i and j in the tangent plane:
```
θᵢⱼ = arccos(dᵢ · dⱼ)
```

For the inner edges to meet, the width at each strut end must be:
```
wᵢ = 2 · L · sin(θᵢⱼ/2) / sin(θᵢⱼ)
```

where L is a characteristic length depending on the triangle geometry.

**Critical insight:** The required width is **different at each vertex** because the angles vary across the dome.

---

## Design Options

### Option 1: Fixed Width + Hub Connectors

- Use uniform strut width across the dome
- Accept that inner edges won't meet perfectly
- Fill gaps with hub connectors (complex 3D shapes)

**Pros:** Simpler struts (all same width)
**Cons:** Complex hub geometry, many unique hub shapes

### Option 2: Tapered Struts + Perfect Inner Triangles

- Compute the exact width needed at each strut end
- Struts are **tapered** (different width at each end)
- Inner faces meet perfectly, forming triangles
- Only need simple "outer caps" to cover external gaps

```
    Narrow end                    Wide end
        │                             │
        ▼                             ▼
    ┌───────┐                   ┌───────────┐
    │       │                   │           │
    │       │───────────────────│           │
    │       │                   │           │
    └───────┘                   └───────────┘
    
    Strut tapers from one vertex to the other
```

**Pros:** Clean inner surface, simpler connectors
**Cons:** Each strut has unique dimensions, more complex cutting

### Option 3: Flat Gusset Plates

- Use uniform cuboid struts
- Connect with flat triangular plates on inner and outer surfaces
- Struts sandwiched between plates, bolted through

```
    ═══════════════════  ← Outer gusset plate
         │   │   │
         │   │   │       ← Struts
         │   │   │
    ═══════════════════  ← Inner gusset plate
```

**Pros:** Simple 2D plates (easy CNC), uniform struts
**Cons:** Visible plates, requires through-bolts

### Option 4: Return to Wedged Struts

- Wedged struts naturally meet at vertices
- No hub needed
- Each strut is a tapered box (CNC-cuttable)

**Pros:** Elegant geometry, no connectors needed
**Cons:** More complex strut shape

---

## Summary

| Approach | Strut Shape | Hub/Connector | Manufacturability |
|----------|-------------|---------------|-------------------|
| Wedged struts | Tapered (radial faces) | None needed | Moderate |
| Cuboid + complex hub | Uniform rectangle | Complex 3D polyhedra | Difficult |
| Tapered cuboid | Tapered rectangle | Simple outer caps | Moderate |
| Cuboid + gussets | Uniform rectangle | Flat plates | Easy |

### The Mathematical Reality

**Cuboid struts with uniform width cannot meet cleanly at geodesic dome vertices.** The geometry fundamentally doesn't allow it because:

1. Strut angles vary across the dome
2. The angles are not coplanar at vertices
3. No simple prismatic hub shape works universally

The solutions are:
- Accept complexity in the hubs (Option 1)
- Move complexity to the struts via tapering (Option 2)
- Use a different connection strategy (Option 3)
- Use wedged struts that naturally meet (Option 4)

---

## Appendix: Angle Analysis for a 3V Honeycomb Dome

In a 3V honeycomb dome, typical angles between struts at vertices:

| Vertex Type | Angle 1 | Angle 2 | Angle 3 | Coplanar? |
|-------------|---------|---------|---------|-----------|
| Type A | 108° | 121° | 121° | No |
| Type B | 118° | 118° | 118° | Nearly |

The variation in angles (108° to 121°) is why no single hub shape works for all vertices.

