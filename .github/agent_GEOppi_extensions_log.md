# GEOppi Extensions Log

## 2026-06-15

- Added a new downstream-building aggregation function in `geoppi/suitable_network_routing.py` to compute cumulative connected buildings for each line segment using a producer-to-consumer graph orientation.
- Exported the new function from `geoppi/__init__.py` so it is directly usable from the package interface.
- Reason for change: support line-segment analysis of connected consumers downstream of producer locations in a radial network context.
