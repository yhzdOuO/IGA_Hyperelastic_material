# IGA cantilever visualization

The module is split into a Qt-independent sampling/evaluation layer and a Qt/OpenGL 3.3 renderer:

- `FieldEvaluator` reuses the solver's arbitrary-point NURBS evaluation and generic hyperelastic material interface, then pushes PK2 stress forward to Cauchy stress.
- `BoundarySampler` samples every nonzero knot span on all six parameter faces and produces consistently outward-wound triangles plus element-boundary line indices.
- `SurfaceMesh` keeps true stress values separate from the draw-only deformation scale.
- `IgaVisualizationWidget` uploads interleaved position/normal/scalar vertices and renders the surface, colorbar, wireframe, and reference configuration.

Run with `--samples-per-span N` (`N >= 2`, default `4`). Controls: left-drag orbit, right-drag pan, wheel zoom, `W` element wireframe, `R` reference configuration, `+/-` deformation scale, `1..8` stress field, `C` palette, and `Esc` close.
