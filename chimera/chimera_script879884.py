try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np

	di=np.array([2652,2654,2656,2658,2662,2665,2666,2667,2668,2669,\
		2670]).astype("uint32")

	ids=np.array([2665,2666,2667,2668]).astype("uint32")

	r=np.array([0.901248,0.901248,0.901248,0.901248]).astype("float")

	clr=np.array([[209,179,139,255],
		[209,179,139,255],
		[209,179,139,255],
		[209,179,139,255]]).astype("uint8")
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
	rc(session,"transparency ~@CD,HD1,HD2,HD3 70 atoms")
	rc(session,"save /Users/albertsmith/Documents/Dynamics/MF_MD_theory/Figures/methyl_rot/fr0_det2.png width 600 height 500 supersample 2 transparentBackground true")
except:
	pass
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script879884.py")
	rc(session,"exit")
