from chimera import runCommand as rc
try:
	rc("open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/1d3z.pdb_5.pdb")
	rc("set bg_color white")
	rc("~ribbon")
	rc("~display")
	rc("display @N,C,CA,O,H,HN")
	rc("set bgTransparency")
	rc("~sel")
except:
	pass