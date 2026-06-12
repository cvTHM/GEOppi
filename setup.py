"""
GEOppi - Conversion of GEO-referenced Piping data into calculable heating networks
"""

from setuptools import setup, find_packages

setup(
    name="GEOppi",
    version="0.2.1",
    description="GEOppi is an open-source software tool designed for different processes in the planning and dimensioning stage of heating networks. It makes use of low-level geo-referenced data for network routings (line objects), heat consumers and heat suppliers (polygon objects) and supports planning processes for new heating networks (e.g. in municipal heat planning) as well as the digital representation of existing networks.",
    packages = find_packages(),
    package_data = {
        'geoppi':['examples/**/*', 'examples/data/exampleNetwork/*', 'QGIS_ModelDesigner/*']
    },
    include_package_data = True,
    install_requires = [
        'numpy',
        'pandas<=2.3.3',
        'geopandas',
        'pandapower',
        'pandapipes<=0.11',
        'rasterio==1.3.10',
        'rasterstats',
        'networkx',
        'libpysal',
    ],
    python_requires="~=3.12"
)
