try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np
	rc(session,"open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/pdbr2000000.pdb_0.pdb")
	rc(session,"~display")
	rc(session,"~ribbon")
	if len(session.models)>1:
		atoms=session.models[1].atoms
		rc(session,"display #1.1")
	else:
		atoms=session.models[0].atoms
		rc(session,"display #1")
	hide=getattr(atoms,"hides")
	hide[:]=1
	hide[di]=0
	setattr(atoms,"hides",hide)
except:
	pass
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script924449.py")
