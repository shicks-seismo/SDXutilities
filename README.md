# SDXutilities
1. Convert multiple SDX files to single QuakeML catalog using ObsPy inventory
Example: 
```
#!/python3
from sdx_utils import sdxtoquakeml
sdxtoquakeml("examples", "test.xml")
```
This will output an xml formatted file named "test.xml"
