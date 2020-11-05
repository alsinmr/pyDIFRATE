from chimera import runCommand as rc
try:
	rc("open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/run1.part0072.xtc_0.pdb")
	rc("set bg_color white")
	rc("~ribbon")
	rc("display ~solvent")
	rc("set bgTransparency")
	rc("~sel")
except:
	pass