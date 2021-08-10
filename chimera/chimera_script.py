try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np

	di=np.array([2308,2310,2312,2314,2318,2321,2322,2323,2324,2325,\
		2326]).astype("uint32")

	ids=np.array([2321,2322,2323,2324]).astype("uint32")

	r=np.array([4.753146,4.753146,4.753146,4.753146]).astype("float")

	clr=np.array([[142,89,77,255],
		[142,89,77,255],
		[142,89,77,255],
		[142,89,77,255]]).astype("uint8")
	mdl=session.open_command.open_data("/Users/albertsmith/Documents/Dynamics/MF_MD_theory/Figures/python/ILE.cxs")[0]
	session.models.add(mdl)
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
	rc(session,"style ball")
	rc(session,"size stickRadius 0.2")
	rc(session,"color all tan")
	r0=getattr(atoms,"radii").copy()
	clr0=getattr(atoms,"colors").copy()
	r0[:]=.8
	r0[ids]=r
	clr0[ids]=clr
	setattr(atoms,"radii",r0)
	setattr(atoms,"colors",clr0)
	rc(session,"")
	rc(session,"save /Users/albertsmith/Documents/Dynamics/MF_MD_theory/Figures/methyl_rot/fr4_det5.png width 600 height 500 supersample 2 transparentBackground true")
except:
	pass
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script727095.py")
	rc(session,"exit")
