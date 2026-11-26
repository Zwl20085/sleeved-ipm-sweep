"""Reference oracle: parse a Maxwell ``.aedt`` project and replay its 2-D geometry history.

Modules: :mod:`.aedt_parse` (file format), :mod:`.aedt_expr` (expression evaluator with
units), :mod:`.aedt_vars` (variable tables), :mod:`.geom2d` (exact 2-D lines/arcs/curves and
Booleans via shapely with exact recovery of every vertex), :mod:`.aedt_geom` (history
replay), :mod:`.model` (loading and evaluating one design). Used only to generate the
reference fixtures in ``tests/fixtures/maxwell_ref``; needs ``shapely``.
"""
