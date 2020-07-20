import os
from chimera import runCommand as rc
try:
	rc("open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/run1.part0002.xtc_1000.pdb")
	rc("set bg_color white")
	rc("~ribbon")
	rc("~display")
	rc("display")
	rc("represent wire")
	rc("linewidth 5")
	rc("color grey")
	rc("defattr /Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/attr_rho183794.txt raiseTool false")
	rc("defattr /Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/attr_radius183794.txt raiseTool false")
	rc("rangecolor rho 0 yellow 0.25 orange 1.0 red")
	rc("sel @/serialNumber=40 or serialNumber=42 or serialNumber=42 or serialNumber=92 or serialNumber=92 or serialNumber=95 or serialNumber=95 or serialNumber=98 or serialNumber=98 or serialNumber=101 or serialNumber=101 or serialNumber=104 or serialNumber=104 or serialNumber=107 or serialNumber=107 or serialNumber=110 or serialNumber=110 or serialNumber=113 or serialNumber=113 or serialNumber=116 or serialNumber=116 or serialNumber=119 or serialNumber=119 or serialNumber=122 or serialNumber=122 or serialNumber=125 or serialNumber=125 or serialNumber=128 or serialNumber=128 or serialNumber=131 or serialNumber=131 or serialNumber=131 or serialNumber=31 or serialNumber=33 or serialNumber=33 or serialNumber=45 or serialNumber=45 or serialNumber=48 or serialNumber=48 or serialNumber=51 or serialNumber=51 or serialNumber=54 or serialNumber=54 or serialNumber=57 or serialNumber=57 or serialNumber=60 or serialNumber=60 or serialNumber=63 or serialNumber=65 or serialNumber=67 or serialNumber=67 or serialNumber=70 or serialNumber=70 or serialNumber=73 or serialNumber=73 or serialNumber=76 or serialNumber=76 or serialNumber=79 or serialNumber=79 or serialNumber=82 or serialNumber=82 or serialNumber=85 or serialNumber=85 or serialNumber=88 or serialNumber=88 or serialNumber=88 or serialNumber=41 or serialNumber=43 or serialNumber=44 or serialNumber=93 or serialNumber=94 or serialNumber=96 or serialNumber=97 or serialNumber=99 or serialNumber=100 or serialNumber=102 or serialNumber=103 or serialNumber=105 or serialNumber=106 or serialNumber=108 or serialNumber=109 or serialNumber=111 or serialNumber=112 or serialNumber=114 or serialNumber=115 or serialNumber=117 or serialNumber=118 or serialNumber=120 or serialNumber=121 or serialNumber=123 or serialNumber=124 or serialNumber=126 or serialNumber=127 or serialNumber=129 or serialNumber=130 or serialNumber=132 or serialNumber=133 or serialNumber=134 or serialNumber=32 or serialNumber=35 or serialNumber=34 or serialNumber=47 or serialNumber=46 or serialNumber=50 or serialNumber=49 or serialNumber=53 or serialNumber=52 or serialNumber=56 or serialNumber=55 or serialNumber=59 or serialNumber=58 or serialNumber=62 or serialNumber=61 or serialNumber=64 or serialNumber=66 or serialNumber=69 or serialNumber=68 or serialNumber=72 or serialNumber=71 or serialNumber=75 or serialNumber=74 or serialNumber=78 or serialNumber=77 or serialNumber=81 or serialNumber=80 or serialNumber=84 or serialNumber=83 or serialNumber=87 or serialNumber=86 or serialNumber=90 or serialNumber=89 or serialNumber=91")
	rc("display sel")
	rc("represent bs sel")
	rc("set bgTransparency")
	rc("~sel")
except:
	print("Error occured in chimera script")
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/attr_radius183794.txt")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/attr_rho183794.txt")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script183794.py")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script183794.pyc")
