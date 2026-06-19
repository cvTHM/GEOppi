# Changelog for GEOppi

## \[v0.2.1]
* \[ADDED] Example script for application of *sum\_attrs\_to\_closest\_supplier*
* \[CHANGED] Function *suitable\_network\_routing* - *sum\_attrs\_to\_closest\_supplier* now takes multiple attributes from attached buildings to be summed up along the shortest path to closest supplier.
* \[CHANGED] Function *geodata\_from\_geometry* for creating pandapipes-compatible geodata from geometry of GIS data supports polygons
* \[CHANGED] Deleted duplicate function definitions in *aux\_functions.py*
* \[CHANGED] Updated setup.py and init file with version number

## \[v0.2.0]

* \[CHANGED] Readme with edited hints for installation
* \[ADDED] Script Examples/ *line\_density\_calculation\_example.ipynb* 
* \[CHANGED] Function *closest\_objects\_to\_points* taking addiditonal argument for returning multiple closest geometries to points within given distance
* \[ADDED] Function *nearest\_points* returning closest points within two geometries
* \[ADDED] Function *get\_adjacent\_lines* to determine closest and resident lines ("Anliegerstraßen")
* \[CHANGED] Function *sum\_attributes\_on\_lines* to optionally take dictionary for matching of polygons and nearest lines instead of determining closest line object

