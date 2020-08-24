import os
from chimera import runCommand as rc
try:
	rc("open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/run1.part0072.xtc_0.pdb")
	rc("set bg_color white")
	rc("~ribbon")
	rc("display")
	rc("set bgTransparency")
	rc("~sel")
	rc("open bild:/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/tensors382714.bild")
	rc("setattr m stickScale .75")
except:
	print("Error occured in chimera script")
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/tensors382714.bild")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script382714.py")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script382714.pyc")
