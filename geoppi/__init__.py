"""
GEOppi - Conversion of GEO-referenced Piping data into calculable heating networks
"""

__version__ = "0.0.1"

# Importing key components
from .auxFunctions import *
from .internal_auxFunctions import *

from .suitable_network_routing import (network_span_bfs, sum_heat_demands_to_closest_supplier, )

from .create_network_topology import (create_ppi_network_from_gdf, )

from .dimension_network_pipes import (update_ppi_results, assign_insulation_type, hydraulicDimensioningNetwork_singleLoadPoint, assign_nominal_widths_manually, )

from .pipe_characteristics import (load_hydraulic_pipe_characteristics, )

