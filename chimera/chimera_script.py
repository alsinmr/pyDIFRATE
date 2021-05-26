try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np

	di=np.array([1354,1356,1358,1359,1360,1361,1362,1363]).astype("uint32")
	rc(session,"open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/HETs_MET_4pw_final_BBfit_1ms.xtc_102001.pdb")
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
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script707742.py")
