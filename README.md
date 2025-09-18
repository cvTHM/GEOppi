# **GEOppi** - Conversion of GEO-referenced Piping data into calculable heating networks
Python library for heating network planning and generation of hydraulic-thermal calculations in [pandapipes](https://github.com/e2nIEE/pandapipes).

![License](https://img.shields.io/badge/license-MIT-green.svg)


**GEOppi** is an open-source software tool designed for different processes in the planning and dimensioning stage of heating networks. It makes use of low-level geo-referenced data for network routings (line objects), heat consumers and heat suppliers (polygon objects) and supports planning processes for new heating networks (e.g. in municipal heat planning) as well as the digital representation of existing networks.

## Overview

**GEOppi** supports planning processes for heating networks in an early and advanced stage. As an open-source tool it is designed to be used freely by planers, public services, energy suppliers and academic researchers / students. 
Three main functionalities ranging from coarse to more granular design phases are implemented which define the possible use cases of the tool:
- Determination of suitable network routings in terms of most profitable line attributes (e.g. line densities / linear heat densities)
- Determination of possible networks expansion based on defined locations and available thermal energy and power budgets of heat supply sites
- Conversion of geo-referenced data of network routings, heat consumers, heat suppliers and periphery into calculable thermal-hydraulic network models using the [pandapipes](https://github.com/e2nIEE/pandapipes) framework


## Features and functionalities
The implemented functionalities rely on input data of possible network routings and spatially distributed heat demands (heat cadastre data). For heating networks, the network course mainly follows street routings, especially in urban settlements. 

**Determination of suitable network routing**  
- Transfer possible network routings as line objects with assigned line attributes to the function
- Freely specify the line attribute which shall be used to determine the segment weights (e.g. line density)
- Specify start points for the search for suitable network routings or leave it to the automatic choice of most profitable line segments
- Specify thresholds for
    - Min. value of line attribute to add segment to the network
    - Min. average value of line attribute of up to three consecutive segments
    - Combinations of max. segment length and min. line attribute value to skip the segment
- Postprocessing options include
    - Merging of touching networks (for multiple start points)
    - Deletion of line segments with attribute below threshold (in descending order of segment length) while preserving hydraulic connectivity
    - Deletion of networks with too few connected consumers

Figure 

<table>
  <tr>
    <td align="center" width="100%">
      <img src="https://raw.githubusercontent.com/cvTHM/GEOppi/refs/heads/main/docs/images/fig_networkRouting_example_v1.png" alt="search spaces" 
width="100%"/><br>
      <sub><b>Figure 1:</b> Fig1</sub>
    </td>
  </tr>
</table>

(The calculation of line attributes is not (yet) part of **GEOppi**. Aggregating attributes of heat consumers to nearest line segments has to be done beforehand.)

**Determination of possible networks expansion**
- 

