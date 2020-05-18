import os
from chimera import runCommand as rc
try:
	rc("open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/1ubi.pdb_0.pdb")
	rc("set bg_color white")
	rc("~ribbon")
	rc("~display")
	rc("display @N,C,CA,O,H,HN")
	rc("represent wire")
	rc("linewidth 5")
	rc("color grey")
	rc("defattr /Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/attr_rho995501.txt raiseTool false")
	rc("defattr /Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/attr_radius995501.txt raiseTool false")
	rc("rangecolor rho 0 yellow 0.25 orange 1.0 red")
	rc("sel :0,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,37,38,39,40,41,42,43,44,45,47,48,49,50,51,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73&@C,O/pdbSegment=A|:1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,36,38,39,40,41,42,43,44,45,46,48,49,50,51,52,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74@N,H,HN,CA/pdbSegment=A")
	rc("display sel")
	rc("represent bs sel")
	rc("set bgTransparency")
	rc("~sel")
except:
	print("Error occured in chimera script")
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/attr_radius995501.txt")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/attr_rho995501.txt")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script995501.py")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script995501.pyc")
