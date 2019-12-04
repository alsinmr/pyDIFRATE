from chimera import runCommand as rc
rc("open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/rmsfit_reduced_1ns.dcd_0.pdb")
rc("set bg_color white")
rc("~ribbon")
rc("~display")
rc("display @N,C,CA,O,H,HN")
rc("set bgTransparency")
rc("~sel")
