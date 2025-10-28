# **GEOppi** - Conversion of GEO-referenced Piping data into calculable heating networks
Python library for heating network planning and generation of hydraulic-thermal calculations in [pandapipes](https://github.com/e2nIEE/pandapipes).

![License](https://img.shields.io/badge/license-MIT-green.svg)


**GEOppi** is an open-source software tool designed for different processes in the planning and dimensioning stage of heating networks. It makes use of low-level geo-referenced data for network routings (line objects), heat consumers and heat suppliers (polygon objects) and supports planning processes for new heating networks (e.g. in municipal heat planning) as well as the digital representation of existing networks.

## Overview

**GEOppi** supports planning processes for heating networks in an early and advanced stage. As an open-source tool it is designed to be used freely by planers, public services, energy suppliers and academic researchers / students. 
Three main functionalities ranging from coarse to more granular design phases are implemented which define the possible use cases of the tool:
- Determination of suitable network routings in terms of most profitable line attributes (e.g. line densities / linear heat densities)
- Determination of possible networks expansion based on defined locations and available thermal energy and power budgets of heat supply sites
- Conversion of geo-referenced data of network routings, heat consumers, heat suppliers and periphery into calculable thermal-hydraulic network models using the [pandapipes](https://github.com/e2nIEE/pandapipes) framework.


<table>
  <tr>
    <td align="center" width="100%">
      <img src="https://raw.githubusercontent.com/cvTHM/GEOppi/refs/heads/main/docs/images/fig_GEOppi_functionalities_Overview.png" alt="search spaces" 
width="100%"/><br>
      <sub><b>Figure 1:</b> Summary of workflows in the planning process for heating networks which can be supported by GEOppi. Top left: Line density values as attributes of street segments (prerequisite, input data); Top right: Suitable network routings on the basis of two differnt heat demand scenarios; Bottom right: Possible network expansion for two different thermal power budgets of connected heat suppliers; Bottom right: Final network topology with dimensioned pipes and distribution of volumetric flow at a design load point </sub>
    </td>
  </tr>
</table>


## Features and functionalities
The implemented functionalities rely on input data of possible network routings and spatially distributed heat demands (heat cadastre data). For heating networks, the network course mainly follows street routings, especially in urban settlements. 

**Determination of suitable network routing**  
Suitable network routings within a complete set of all possible network routings (e.g. street system) are found by scanning the given system with a modified depth-first-search (DFS) algorithm. It is applied to the [networkx](https://github.com/networkx/networkx) multigraph representation of the network. The process comes to an end if none of the reachable, i.e. hydraulically connected to the spanned network, branches features line attributes meeting the defined thresholds.
- Transfer all possible network routings as line objects with assigned line attributes to the function (e.g. street system)
- Freely specify the line attribute which shall be used to determine the segment weights (e.g. line density)
- Specify start points for the search for suitable network routings or leave it to the automatic choice of most profitable line segments
- Specify thresholds for
    - Min. value of line attribute to add segment to the network
    - **Special feature**: Min. average value of line attribute of up to three consecutive segments
    - Combinations of max. segment length and min. line attribute value to skip the segment
- Postprocessing options include
    - Merging of touching networks (for multiple start points)
    - Deletion of line segments with attribute below threshold (in descending order of segment length) while preserving hydraulic connectivity
    - Deletion of networks with too few connected consumers

<table>
  <tr>
    <td align="center" width="100%">
      <img src="https://raw.githubusercontent.com/cvTHM/GEOppi/refs/heads/main/docs/images/fig_networkRouting_example_v1.png" alt="search spaces" 
width="30%"/><br>
      <sub><b>Figure 2:</b> Exemplified process for advanced look-up in three consecutive line segments to assess the average line density in the currentlyx scanned branch. </sub>
    </td>
  </tr>
</table>

(The calculation of line attributes is not yet part of **GEOppi**. Aggregating attributes of heat consumers to nearest line segments has to be done beforehand.)

**Determination of possible networks expansion**  
A maximum expansion of a heating network starting at specified heat supply sites is computed by successively subtracting thermal demands of scanned line segments (annual energy demands and thermal power demands at a design load point) from the available energy and power budgets of the heat supply site. Thermal losses in the network and simultaneity factors for the summed thermal power demanded can be considered as constant or dynamically mutable factors. Distances between network branch ends and heat supply sites are kept short with the use of a modified breadth-first-search (BFS) algorithm. It is applied to the [networkx](https://github.com/networkx/networkx) multigraph representation of the network.
- Transfer all possible network routings as line objects with assigned line attributes to the function (e.g. street system)
- Specify locations of heat supply sites as polygons and their attributes (available thermal energy and power budgets)
- **Special feature**: Provide functional relationships for (a) relative thermal losses in the network in dependence of average line density and (b) for simultaneity factors in dependence of the number of connected buildings/heat consumers
- Postprocessing options include
  - Merging of touching networks (for multiple start points)
  - Deletion of line segments with attribute below threshold (in descending order of segment length) while preserving hydraulic connectivity
  - Deletion of networks with too few connected consumers


**Conversion of geo-referenced data into calculable network model**  
Geo-referenced line data for network routing and polygon data for heat consumers/buildings and heat supply sites (+ point objects for valves as network periphery) are converted into calculable network models (hydraulically and thermally with supply line and return line) in the [pandapipes](https://github.com/e2nIEE/pandapipes) framework.  
The network topology is modelled as a node-edge-representation. Heat consumers and heat supply sites are modelled using standard components from the [pandapipes component compendium](https://pandapipes.readthedocs.io/en/latest/components.html) for versions >= 0.10.0.  
An **automatic pipe dimensioning routine** is selectable - pipe data for commercially available plastic jacket pipes and simple PE-HD pipes are implemented.

<table>
  <tr>
    <td align="center" width="100%">
      <img src="https://raw.githubusercontent.com/cvTHM/GEOppi/refs/heads/main/docs/images/fig_GEOppi_exampleNetwork.png" alt="search spaces" 
width="100%"/><br>
      <sub><b>Figure 3:</b> Chronological steps of heating network construction from GIS data in **GEOppi** </sub>
    </td>
  </tr>
</table>

- Load polygon objects of heat consumers, preferrably with attributes for heat demand
- Load polygon objects of heat supply sites, preferrably with characteristics for type (base load producer or peak load producer) and thermal power
- (Optional) Load point objects of valves as network periphery
- (Optional) Provide rasterized data of digital elevation model for geodetic heights at nodes
- Load line objects of the final network topology (distribution lines and house connection lines)
  - Supports completely new heating networks (no attributes on the line objects) and existing heating network topologies (specify attributes which shall be maintained)
  - Only the supply line layer is needed for this step!
  - Note: The automatic creation of house connection lines for entirely new heating networks is not yet part of **GEOppi**. For this reason, a dedicated function for the [QGIS model designer](https://docs.qgis.org/3.40/en/docs/user_manual/processing/modeler.html) is implemented in the [QGIS_ModelDesigner directory](https://github.com/cvTHM/GEOppi/blob/main/geoppi/QGIS_ModelDesigner). It creates straight house connection lines from each selected building in the masked investigation area to the closest distribution line segment and a robust separation of the distribution line segment at the connection point.
- For the case of entirely new heating networks:
  - Specify type of piping (plastic jacket pipes (+ stage of insulation) or simple PE-HD pipes (uninsulated))
  - **Special feature**: For automatic pipe dimensioning routine, specify:
    - A design load point defining thermal powers at each heat consumer and the design temperature spread in the network
    - A target specific pressure loss
    - The delivered themal power of each base load producer at the design load point
    - Sets of nominal widths which are allowed for distribution lines and for house connection lines
  - Postprocessing options include
    - Assignment of nominal widths and corresponding hydraulic and thermal pipe characteristics to a set of selected pipes simultaneously
    - Detection of cycles in a meshed network
    - **Special feature**: Implementation of controllers in final networks (supply and return line)
      - Control of fixed return temperature at selectable heat consumers
      - Control of min. absolute pressure in the network
      - Control of min. differential pressure at worst-supplied heat consumers
      - Control of cascaded operation of heat suppliers (specified merit order)

## Technical requirements
**GEOppi** requires the following essential modules in the used Python environment

- [pandapipes](https://github.com/e2nIEE/pandapipes) version >= 0.10.0
- [pandapower](https://github.com/e2nIEE/pandapower)
- [NumPy](https://github.com/numpy/numpy)
- [Pandas](https://github.com/pandas-dev/pandas)
- [GeoPandas](https://github.com/geopandas/geopandas)
- [Rasterio](https://github.com/rasterio/rasterio)
- [networkx](https://github.com/networkx/networkx)

**Used input data types**
- Geo-referenced data are taken as inputs and given as outputs in GPKG format (applies to all geometry objects)

## Installation
You can install **GEOppi** by cloning a branch or copying the compressed branch directory to your local machine and running the enclosed setup file, e.g. in a shell after directing to the target folder:

```
python -m pip install .
```

This command installs **GEOppi** and its main dependencies. 
Please note that problems may still occur when trying to install [GeoPandas](https://github.com/geopandas/geopandas) or [Rasterio](https://github.com/rasterio/rasterio) via pip. You can manually install them if the installation fails during the process for **GEOppi**.

## Quick start
Please refer to the template scripts [networkExtension_example.py](https://github.com/cvTHM/GEOppi/blob/main/geoppi/examples/networkExtension_example.py) and [networkModelling_example.py](https://github.com/cvTHM/GEOppi/blob/main/geoppi/examples/networkModelling_example.py) for a quick start.

## Documentation
A central documentation has not been started yet. Exemplary applications of the functionalities for determination of possible network expansions and the conversion of geo-referenced data into a calculable network model are given in the files [networkExtension_example.py](https://github.com/cvTHM/GEOppi/blob/main/geoppi/examples/networkExtension_template.py) and [networkModelling_example.py](https://github.com/cvTHM/GEOppi/blob/main/geoppi/examples/networkModelling_template.py).

## Contributing
All contributions to the tool **GEOppi** are warmly welcomed. Refer to the [contribution guidelines](https://github.com/cvTHM/GEOppi/blob/main/CONTRIBUTION_guidelines.md) for more information.

## License
The software tool is licensed under the [MIT License](https://github.com/cvTHM/GEOppi/blob/main/LICENSE).

## Citation
**GEOppi** and its functionalities are first introduced to public in a conference article.
When using **GEOppi** in your research, please cite it as follows:

Völzel, C., Schug, N., Textor, M., Lechner, S. (2025): Constructing heating networks from GIS data: A framework for editing and converting line topologies into calculable grids. Konferenzbeitrag, NEIS 2025 Conference on Sustainable Energy Supply and Energy Storage, Hamburg, September 2025

## Contact
Please refer to the issues or the discussion section on the GitHub repository for questions and feedback.
