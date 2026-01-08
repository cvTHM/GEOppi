"""
GEOppi - Conversion of GEO-referenced Piping data into calculable heating networks
"""

from setuptools import setup, find_packages

setup(
    name="GEOppi",
    version="0.1.0",
    description="GEOppi is an open-source software tool designed for different processes in the planning and dimensioning stage of heating networks. It makes use of low-level geo-referenced data for network routings (line objects), heat consumers and heat suppliers (polygon objects) and supports planning processes for new heating networks (e.g. in municipal heat planning) as well as the digital representation of existing networks.",
    packages = find_packages(),
    package_data = {
        'geoppi':['examples/**/*', 'examples/data/exampleNetwork/*', 'QGIS_ModelDesigner/*']
    },
    include_package_data = True,
    install_requires = [
        'numpy',
        'pandas',
        'geopandas',
        'pandapower==2.14.11',
        'pandapipes>=0.10.0',
        'rasterio',
        'networkx',
    ],
    python_requires=">=3.10"
)
