import os
import numpy as np
from matplotlib.tri import Triangulation
from chimerax.core.models import Surface
from chimerax.surface import calculate_vertex_normals,combine_geometry_vntc
from chimerax.core.commands import run as rc

def surf3D(session,x,y,z,colors=None,triangles=None):
    """
    Creates a 3D surface plot in chimera. Input are the x,y, and z coordinates
    (x and y can be the same size as z, or may be x and y axes for the z plot, 
    in which case x.size==z.shape[0], y.size==z.shape[1].)
    
    Colors may also be specified for each vertex (should then have 3xN or 4xN
    points, where N=z.size), or a single color may be specified for the whole 
    plot. 4xN allows specification of the opacity. Us 0-255 RGB specification.
    """
    
    if colors is None:
        colors=np.array([[210,180,140,255]]).T.repeat(z.size,axis=1)
    else:
        colors=np.array(colors)
        if colors.size==3:colors=np.concatenate((colors,[255]))
        if colors.size==4:colors=np.atleast_2d(colors).T.repeat(z.size,axis=1)
        if colors.ndim>2:colors=np.reshape(colors,[colors.shape[0],colors.shape[1]*colors.shape[2]])
    
    if not(x.size==z.size and y.size==z.size):
        x,y=[q.reshape(x.size*y.size) for q in np.meshgrid(x,y)]
    if triangles is None:
        triangles=Triangulation(x,y).triangles
    z=z.reshape(z.size)
    
    xyz=np.ascontiguousarray(np.array([x,y,z]).T,np.float32)       #ascontiguousarray forces a transpose in memory- not just editing the stride
    colors = np.array(colors, np.uint8).T
    tri = np.array(triangles, np.int32)
    
    norm_vecs=calculate_vertex_normals(xyz,tri)
    
    s=Surface('surface',session)
    s.set_geometry(xyz,norm_vecs,tri)
    s.vertex_colors=colors
    session.models.add([s])


try:

	x0=np.array([-35.000000,-34.861809,-34.723618,-34.585427,-34.447236,-34.309045,-34.170854,-34.032663,-33.894472,-33.756281,\
		-33.618090,-33.479899,-33.341709,-33.203518,-33.065327,-32.927136,-32.788945,-32.650754,-32.512563,-32.374372,\
		-32.236181,-32.097990,-31.959799,-31.821608,-31.683417,-31.545226,-31.407035,-31.268844,-31.130653,-30.992462,\
		-30.854271,-30.716080,-30.577889,-30.439698,-30.301508,-30.163317,-30.025126,-29.886935,-29.748744,-29.610553,\
		-29.472362,-29.334171,-29.195980,-29.057789,-28.919598,-28.781407,-28.643216,-28.505025,-28.366834,-28.228643,\
		-28.090452,-27.952261,-27.814070,-27.675879,-27.537688,-27.399497,-27.261307,-27.123116,-26.984925,-26.846734,\
		-26.708543,-26.570352,-26.432161,-26.293970,-26.155779,-26.017588,-25.879397,-25.741206,-25.603015,-25.464824,\
		-25.326633,-25.188442,-25.050251,-24.912060,-24.773869,-24.635678,-24.497487,-24.359296,-24.221106,-24.082915,\
		-23.944724,-23.806533,-23.668342,-23.530151,-23.391960,-23.253769,-23.115578,-22.977387,-22.839196,-22.701005,\
		-22.562814,-22.424623,-22.286432,-22.148241,-22.010050,-21.871859,-21.733668,-21.595477,-21.457286,-21.319095,\
		-21.180905,-21.042714,-20.904523,-20.766332,-20.628141,-20.489950,-20.351759,-20.213568,-20.075377,-19.937186,\
		-19.798995,-19.660804,-19.522613,-19.384422,-19.246231,-19.108040,-18.969849,-18.831658,-18.693467,-18.555276,\
		-18.417085,-18.278894,-18.140704,-18.002513,-17.864322,-17.726131,-17.587940,-17.449749,-17.311558,-17.173367,\
		-17.035176,-16.896985,-16.758794,-16.620603,-16.482412,-16.344221,-16.206030,-16.067839,-15.929648,-15.791457,\
		-15.653266,-15.515075,-15.376884,-15.238693,-15.100503])

	