# SDXutilities
1. Convert multiple SDX files to single QuakeML catalog using ObsPy inventory
Example: 
```
from sdx_utils import sdxtoquakeml
sdxtoquakeml("examples", "tmp.xml")
```
This will output an xml formatted file named "tmp.xml"
