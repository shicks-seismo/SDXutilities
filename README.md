# SDXutilities
1. sdxtoquakeml: Convert multiple SDX files to single QuakeML catalog using ObsPy inventory
Example: 
```
#!/python3
from sdx_utils import sdxtoquakeml
sdxtoquakeml("examples", "test.xml")
```
This will output an xml formatted file named "test.xml"

2. quakemltosdx: Convert single QuakeML catalog to SDX file for each event
Example:
```
#!/python3
from sdx_utils import quakemltosdx
quakemltosdx("test.xml",
             "stations.dat",
             ["sdx_crust1.txt", 1.73, 2, 10, 20, 20, 500])
```
stations.dat is a text file containing station_code, latitude(dd), longitude(dd) and elevation in m.
