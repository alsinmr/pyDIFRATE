#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 13:17:18 2020

@author: albertsmith
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:05:03 2020

@author: albertsmith

Collection of functions for plotting residual tensors onto molecules in Chimera
"""


import os
curdir=os.getcwd()
import numpy as np
from scipy.spatial import Delaunay
from mpl_toolkits import mplot3d
#from mayavi.mlab import triangular_plot
import matplotlib.pyplot as plt
import matplotlib.colors as colors
os.chdir('../Struct')
from vf_tools import d2,Spher2pars,norm,getFrame,Rspher,sc2angles
os.chdir('../chimera')
from chimera_funs import chimera_path,run_command
from shutil import copyfile
os.chdir(curdir)


def get_path(filename=None):
    """This function opens a file for writing, putting it in the same folder as 
    the chimera_funs.py script"""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    if filename is None:
        return dir_path
    else:
        full_path=os.path.join(dir_path,filename)
        return full_path
    
    
def tensor2xyz(delta,eta,alpha,beta,gamma,sc=2.18,q=8):
    """
    Calculates x,y, and z coordinates for a rank-2 tensor having delta and eta, 
    and euler angles alpha, beta, and gamma. The scaling is set so that an 
    unaveraged tensor spans the full length of an H–C bond. Can be adjusted for
    other bond lengths.
    """
    
    sc=np.sqrt(2/3)*sc
    
    a,b=sphere(q)
    
    A=[-1/2*delta*eta,0,np.sqrt(3/2)*delta,0,-1/2*delta*eta]   #Components in PAS
    
    #0 component after rotation by a and b
    A0=np.array([A[mp+2]*d2(b,m=0,mp=mp)*np.exp(1j*mp*a) for mp in range(-2,3)]).sum(axis=0).real
    
    #Coordinates before rotation by alpha, beta, gamma
    x0=np.cos(a)*np.sin(b)*np.abs(A0)*sc/2
    y0=np.sin(a)*np.sin(b)*np.abs(A0)*sc/2
    z0=np.cos(b)*np.abs(A0)*sc/2

    
    #Rotate by alpha
    x1,y1,z1=x0*np.cos(alpha)+y0*np.sin(alpha),-x0*np.sin(alpha)+y0*np.cos(alpha),z0
    #Rotate by beta
    x2,y2,z2=x1*np.cos(beta)-z1*np.sin(beta),y1,np.sin(beta)*x1+np.cos(beta)*z1
    #Rotate by gamma
    x,y,z=x2*np.cos(gamma)+y2*np.sin(gamma),-x2*np.sin(gamma)+y2*np.cos(gamma),z2

    """
    Here, we create a group index. The plots come out nicer if we separate positive
    and negative parts (A0) and separate different lobes of the tensor
    """
    
#    index=list()
#    i=np.argwhere(A0>=0).squeeze()
#    index.append(i)
#    index.append(i+n)
#    
#    if eta>0.15:
#        i=np.argwhere(np.logical_and(A0<=0,np.logical_or(a<=np.pi/2,a>3*np.pi/2))).squeeze()
#        index.append(i)
#        i=np.argwhere(np.logical_and(A0<=0,np.logical_or(a>np.pi/2,a<=3*np.pi/2))).squeeze()
#        index.append(i)
#    else:
#        i=np.argwhere(A0<=0).squeeze()
#        index.append(np.concatenate((i,i+n)))
    
#    index=[np.argwhere(A0>=0).squeeze(),np.argwhere(A0<=0).squeeze()]

    return x,y,z,A0

def sphere(q=8):
    """
    Returns a set of Euler angles (alpha,beta) over a sphere. Note- this is not
    a valid powder average for simulations– the first alpha value is redundant
    with the last alpha value for each beta (that is, alpha[-1]=alpha[0]+2*pi,
    for each beta). This allows us to use a triangulation algorithm for 3D 
    plotting but would distort simulation results.
    
    alpha,beta = sphere(q=8)
    """
    NP=np.ceil(q*np.sin(np.arange(0.5,q+0.5)*np.pi/(2*q))).astype(int)
    NP41=4*NP+1
    
    nF0=NP.sum()
    nF = NP41.sum()
    
    alpha=np.zeros(nF)
    beta=np.zeros(nF)
    
    theta_j = 0;
    count = 0;

    for j in range(q):
        dtheta = np.arccos(np.cos(theta_j) - NP[j]/nF0 ) - theta_j
        beta[np.arange(count,count+NP41[j])] = theta_j + dtheta/2
        dphi = np.pi/(2*NP[j]);
        alpha[count:count+4*NP[j]] = np.linspace(0.5*dphi,(4*NP[j] - 0.5)*dphi,4*NP[j]);
        alpha[count+4*NP[j]]=alpha[count+4*NP[j]-1]-2*np.pi
        count = count + NP41[j];
        theta_j += dtheta;
    
    alpha=np.concatenate(([0],alpha,alpha[::-1],[0]))
    beta=np.concatenate(([0],beta,np.pi-beta[::-1],[np.pi]))
    
    return alpha,beta
            

def draw_tensorsX(data,tensors='A_m0_finF',index=None,tstep=0,sc=2.09,fileout=None,\
                  scene=None,pos=None,png_opts=None,chimera_cmds=None,marker=None):
    """
    Plots averaged tensors onto a molecule in chimera. A molecule object is required,
    where sel1/sel2 determine which bonds the molecules are to be plotted onto.
    The tensors themselves should also be provided, where they are defined in
    the frame of the bond. These may be provided by either giving a data object
    as the first argument, which then contains both the molecule objection and
    the tensors (the tensors should be in data.vars, and then one provides a
    string 'tensors' which gives the appropriate key to find the tensors). 
    Alternatively, one may provide the molecule object as the first argument,
    in which case 'tensors' should be a 5xN numpy array with the tensors given 
    directly.
    
    Note, the molecule object must still contain the tensor frame function 
    (mol._vft is not None). This should return the orientation of the bond. 
    mol.sel1 and mol.sel2 should define the atoms yielding the bond position.
    
    Optional arguments:
        sc:     Deterines the length of tensors. If provided, this value specifies
                the length of a tensor with delta=1,eta=0. If not provide or None,
                the largest tensor is scaled to extend 2x the mean bond length,
                which is determined from the atom positions in data.sens.mol.sel1
                and data.sens.mol.sel2 (if pos is provided, then bond lengths 
                are instead assumed to be 1.09)
        q:      Quality of the tensor plot. Default is 10 (smoother tensors for
                higher values)
        index:  Index of the tensors to be plotted. Default plots all tensors
        pos:    Position of the tensors to plotted. Default position taken from
                mean position between mol.sel1 and mol.sel2
        tensors:Specify which tensor to plot or provide the tensors directly. 
                Options are "D2inf" (D2 evaluated at infinite time), "avg_tensor"
                which is the averaged rank-2 tensor, or one may provide a dict
                specifying the tensor directly (with delta (N,),eta(N,),euler(3,N)).
                Note, if tensors is provided directly, one may replace the data
                object with a molecule object
        fileout:Save the resulting figure to a png
        scene:  Load a scene saved previously in chimera
        png_opts:A string passed to chimera while saving the image
                (Command: copy file {fileout} png {png_opts})
        chimera_cmds:List of valid chimera commands (executed within chimera)
        marker: Color in one or more tensors with a different color (green/yellow)
                instead of red/blue). List of indices or logical index
                
    draw_tensorsX(data,sc=None,tensors='A_m0_finF',index=None,fileout=None,scene=None,pos=None,png_opts=None,chimera_cmds=None)
    
    """
    
    
    "Setup"
    if not isinstance(tensors,str) and hasattr(data,'pdb'):
        mol=data
    elif isinstance(tensors,str):
        if tensors in data.vars:
            tensors=data.vars[tensors]
        mol=data.sens.molecule
            
    if index is None:
        index=np.arange(tensors.shape[1],dtype=int)
    else:
        index=np.array(index,dtype=bool) if (np.size(index)==tensors.shape[1] and np.max(index)<2) else np.array(index,dtype=int)
        tensors=tensors[:,index]

    

    mol.mda_object.trajectory[tstep]   #Go to the requested time step

    "Get the current positions of the bonds"
    if pos is None:
        pos=(mol.sel1.positions[index]+mol.sel2.positions[index]).T/2
    
    "Get the current orientations of the bonds, rotate tensors accordingly"
    vZ,vXZ=mol._vft()
    vZ=norm(vZ)
    scF=getFrame(vZ[:,index],vXZ[:,index])
    tensors=Rspher(tensors,*scF)
    
    "Convert tensors into parameters"
    tensors=Spher2pars(tensors)
    
    "Make sure we have the pdb or a path to the scene was given"
    if mol.pdb is not None:
        a=np.char.find(mol.pdb[::-1],'.')
        b=np.char.find(mol.pdb[::-1],'_')
        ts=np.array(mol.pdb[-b:-a-1],dtype=int)
        if tstep!=ts and scene is None:
            mol.pdb=None
            
    sel=None
    if scene is not None:
        pdb=None
#    elif mol.pdb is not None:
#        pdb=mol.pdb
    else:
        try:
            if np.unique(mol.sel1.segids).size==1 and np.unique(mol.sel1.resids).size==1:
                select='resid {0} and segid {1}'.format(mol.sel1[0].resid,mol.sel1[0].segid)
            elif np.unique(mol.sel1.segids).size==1:
                select='segid {1}'.format(mol.sel1[0].segid)
            elif np.unique(mol.sel1.resids).size==1:
                select='resid {0}'.format(mol.sel1[0].resid)
            sel=mol.mda_object.select_atoms(select)
            mol.mda_object.trajectory[tstep]
            pdb=get_path('pdb{0}.pdb'.format(tstep))
            mp=sel.positions.mean(0)
            sel.positions=sel.positions-mp
            pos=(pos.T-mp).T
            sel.write(pdb)
#            mol.MDA2pdb(tstep=tstep,select=select)
#            pdb=mol.pdb
        except:
            print('Failed to created pdb for drawing')
            return
    

    delta=tensors[0]
    eta=tensors[1]
    euler=tensors[2:]
        
    theta_steps=100
    phi_steps=50
    positive_color=(255,100,100,255)
    negative_color=(100,100,255,255)
    
    run_chimeraX(delta,eta,euler,pos,theta_steps,phi_steps,positive_color,\
                 negative_color,fileout,png_opts,chimera_cmds,sc=sc,pdb=pdb,scene=scene,marker=marker)
   
def only_tensorsX(tensors=None,delta=None,eta=None,euler=None,index=None,fileout=None,pos=None,png_opts=None,chimera_cmds=None):
    """
    Draws tensors in chimeraX, using python options.
    
    only_tensorsX(delta,eta=None,euler=None,fileout=None,pos=None,png_opts=None)
    
    If multiple tensors are given, this will automatically space them out. One
    may also provide the position
    
    chimera_cmds is a list of strings that are applied as commands in chimera
    
    fileout will save the chimera window as a png. Program will automatically close
    after saving! (intended for scripting multiple images, hence we don't want
    lots of chimera instances staying open)
    
    png_opts is a string appended to the 'save' command if fileout is used
    
    """
    
    
    
    delta=np.atleast_1d(delta)
    n=delta.size
    eta=np.zeros(n) if eta is None else np.atleast_1d(eta)
    euler=np.zeros([3,n]) if euler is None else np.atleast_2d(euler)
    if euler.shape[0]==1:euler=euler.T
    
    "Apply index"
    if index is None:index=np.ones(n,dtype=bool)
    index=np.ones(n,dtype=bool) if index is None else np.atleast_1d(index)
    delta=delta[index]
    eta=eta[index]
    euler=euler[:,index]
    
    
    theta_steps=100
    phi_steps=50
    positive_color=(255,100,100,255)
    negative_color=(100,100,255,255)
    run_chimeraX(delta,eta,euler,pos,theta_steps,phi_steps,positive_color,\
                 negative_color,fileout,png_opts,chimera_cmds)

def run_chimeraX(delta,eta,euler,pos,theta_steps,phi_steps,
                 positive_color,negative_color,fileout,png_opts,
                 chimera_cmds,sc=1,pdb=None,scene=None,marker=None):
    n=delta.size
    rand_index=np.random.randint(1e6)
    full_path=get_path('chimera_script{0:06d}.py'.format(rand_index))
    
    tensor_file=get_path('tensors{0:06d}.txt'.format(rand_index))
    
    "Get the positions of each tensor"
    if pos is None:
        step=np.abs(delta).max()*1.1
        pos=np.zeros([3,n])
        pos[0]=np.linspace(0,n*step,n,endpoint=False)
    else:
        pos=np.atleast_2d(pos)
        if pos.shape[0]==1:pos=pos.T
        
    "Write tensors to file"
    write_tensor(tensor_file,delta,eta,euler,pos,marker)
    

    
    with open(full_path,'w') as f:
        py_line(f,'import os')
        py_line(f,'import numpy as np')
        py_line(f,run_command(version='X'))
        
        copy_funs(f)    #Copy required functions into chimeraX script
                
        
        py_line(f,'try:')
        if scene is not None:
            WrCC(f,'open '+scene,1)
        elif pdb is not None:
            WrCC(f,'open '+pdb,1)
        py_line(f,'load_surface(session,"{0}",sc={1},theta_steps={2},phi_steps={3},positive_color={4},negative_color={5})'\
                    .format(tensor_file,sc,theta_steps,phi_steps,positive_color,negative_color),1)
        if chimera_cmds is not None:
            if isinstance(chimera_cmds,str):chimera_cmds=[chimera_cmds]
            for cmd in chimera_cmds:
                WrCC(f,cmd,1)
        if fileout is not None:
            if fileout[-4:]!='.png':fileout=fileout+'.png'
            if png_opts is None:png_opts=''
            WrCC(f,"save " +fileout+' '+png_opts,1)
            
        py_line(f,'except:')
        py_line(f,'print("Error in chimera script")',1)
        py_line(f,'finally:')
        py_line(f,'os.remove("{0}")'.format(full_path),1)
        py_line(f,'os.remove("{0}")'.format(tensor_file),1)
        if pdb is not None:py_line(f,'os.remove("{0}")'.format(pdb))
        if fileout is not None: #Exit if a file is saved
            WrCC(f,'exit',1)
        
            
    copyfile(full_path,full_path[:-9]+'.py')
    copyfile(tensor_file,tensor_file[:-10]+'.txt')
    
    os.spawnl(os.P_NOWAIT,chimera_path(version='X'),chimera_path(version='X'),full_path)

def copy_funs(f):
    """
    Copys all functions in THIS file below the comment "Files used inside ChimeraX"
    
    Input is the file handle, f, to which the pythons functions should be copied
    
    copy_funs(f)
    """
    
    with open(get_path('vis3D.py'),'r') as funs:
        start_copy=False
        for line in funs:
            if start_copy:
                f.write(line)
            else:
                if len(line)>=30 and line[:30]=="#%% Files used inside ChimeraX":
                    start_copy=True
     
def py_line(f,text,nt=0):
    """
    Prints a line to a file for reading as python code. Inserts the newline and
    also leading tabs (if nt specified)
    
    python_line(f,text,nt=0)
    """
    
    for _ in range(nt):
        f.write('\t')
    f.write(text)
    f.write('\n')

def write_tensor(filename,delta,eta=None,euler=None,pos=None,marker=None):
    """
    Writes out a tab-separated file with delta, eta, alpha, beta, gamma, and
    x,y,z for tensors. For reading within ChimeraX
    
    write_tensor(filename,delta,eta=None,euler=None,pos=None,marker=None)
    """
    
    delta=np.array(delta)
    n=delta.size
    
    #Defaults, make sure all numpy arrays
    eta=np.zeros(n) if eta is None else np.array(eta)
    euler=np.zeros([3,n]) if euler is None else np.array(euler)
    pos=np.zeros([3,n]) if pos is None else np.array(pos)
    if marker is None:
        marker=np.zeros(n)
    else:
        if not(hasattr(marker,'__len__')):marker=[marker]
        if len(marker)<len(eta) or np.max(marker)>1:
            m1=marker
            marker=np.zeros(n)
            marker[np.array(m1,dtype=int)]=1
    
    if len(euler)==3:
        alpha,beta,gamma=euler
    else:
        alpha,beta,gamma=sc2angles(*euler)
    X,Y,Z=pos
    
    with open(filename,'w') as f:
        for vals in zip(delta,eta,alpha,beta,gamma,X,Y,Z,marker):
            for v in vals[:-1]:f.write('{0:16.8}\t'.format(v))
            f.write('{0:d}\t'.format(int(vals[-1])))
            f.write('\n')


def WrCC(f,command,nt=1):
    "Function to print chimera commands correctly"
    for _ in range(nt):
        f.write('\t')
    f.write('rc(session,"{0}")\n'.format(command))
#%% Files used inside ChimeraX (don't edit this comment!!..it will break the code)
"""
Everything after these lines is printed into the chimeraX script, so don't add
anything below that you don't need in chimeraX
"""
def sphere_triangles(theta_steps=100,phi_steps=50):
    """
    Creates arrays of theta and phi angles for plotting spherical tensors in ChimeraX.
    Also returns the corresponding triangles for creating the surfaces
    """
    
    theta=np.linspace(0,2*np.pi,theta_steps,endpoint=False).repeat(phi_steps)
    phi=np.repeat([np.linspace(0,np.pi,phi_steps,endpoint=True)],theta_steps,axis=0).reshape(theta_steps*phi_steps)
    
    triangles = []
    for t in range(theta_steps):
        for p in range(phi_steps-1):
            i = t*phi_steps + p
            t1 = (t+1)%theta_steps
            i1 = t1*phi_steps + p
            triangles.append((i,i+1,i1+1))
            triangles.append((i,i1+1,i1))
    
    return theta,phi,triangles
    
def spherical_surface(delta,eta=None,euler=None,pos=None,sc=2.09,
                      theta_steps = 100,
                      phi_steps = 50,
                      positive_color = (255,100,100,255), # red, green, blue, alpha, 0-255 
                      negative_color = (100,100,255,255)):
    """
    Function for generating a surface in ChimeraX. delta, eta, and euler angles
    should be provided, as well positions for each tensor (length of all arrays
    should be the same, that is (N,), (N,), (3,N), (3,N) respectively.
    
    Returns arrays with the vertices positions (Nx3), the triangles definitions
    (list of index triples, Nx3), and a list of colors (Nx4)
    
    xyz,tri,colors=spherical_surface(delta,eta=None,euler=None,pos=None,
                                     theta_steps=100,phi_steps=50,
                                     positive_color=(255,100,100,255),
                                     negative_color=(100,100,255,255))
    """
    # Compute vertices and vertex colors
    a,b,triangles=sphere_triangles(theta_steps,phi_steps)
    
    if euler is None:euler=[0,0,0]
    if pos is None:pos=[0,0,0]
    if eta is None:eta=0
    
    # Compute r for each set of angles
    sc=np.sqrt(2/3)*sc
    
    A=[-1/2*delta*eta,0,np.sqrt(3/2)*delta,0,-1/2*delta*eta]   #Components in PAS
    
    #0 component after rotation by a and b
    A0=np.array([A[mp+2]*d2(b,m=0,mp=mp)*np.exp(1j*mp*a) for mp in range(-2,3)]).sum(axis=0).real
    
    #Coordinates before rotation by alpha, beta, gamma
    x0=np.cos(a)*np.sin(b)*np.abs(A0)*sc/2
    y0=np.sin(a)*np.sin(b)*np.abs(A0)*sc/2
    z0=np.cos(b)*np.abs(A0)*sc/2

    alpha,beta,gamma=euler
    #Rotate by alpha
    x1,y1,z1=x0*np.cos(alpha)+y0*np.sin(alpha),-x0*np.sin(alpha)+y0*np.cos(alpha),z0
    #Rotate by beta
    x2,y2,z2=x1*np.cos(beta)-z1*np.sin(beta),y1,np.sin(beta)*x1+np.cos(beta)*z1
    #Rotate by gamma
    x,y,z=x2*np.cos(gamma)+y2*np.sin(gamma),-x2*np.sin(gamma)+y2*np.cos(gamma),z2

    x=x+pos[0]
    y=y+pos[1]
    z=z+pos[2]
    
#    xyz=[[x0,y0,z0] for x0,y0,z0 in zip(x,y,z)]
    #Determine colors
    colors=np.zeros([A0.size,4],np.uint8)
    colors[A0>=0]=positive_color
    colors[A0<0]=negative_color
    

    # Create numpy arrays
#    xyz = np.array(xyz, np.float32)
    xyz=np.ascontiguousarray(np.array([x,y,z]).T,np.float32)       #ascontiguousarray forces a transpose in memory- not just editing the stride
    colors = np.array(colors, np.uint8)
    tri = np.array(triangles, np.int32)

    return xyz,tri,colors
 

def load_tensor(filename):
    """
    Reads in a tab-separated file with delta, eta, alpha,beta, gamma, and x,y,z
    for a set of tensors. 
    
    delta,eta,euler,pos=load_tensor(filename)
    """
    delta=list()
    eta=list()
    alpha=list()
    beta=list()
    gamma=list()
    x=list()
    y=list()
    z=list()
    marker=list()
    with open(filename,'r') as f:
        for line in f:
            out=line.strip().split('\t')
            out=[np.array(o,float) for o in out]
            delta.append(out[0])
            eta.append(out[1])
            alpha.append(out[2])
            beta.append(out[3])
            gamma.append(out[4])
            x.append(out[5])
            y.append(out[6])
            z.append(out[7])
            marker.append(out[8])

    delta=np.array(delta)
    eta=np.array(eta)
    euler=np.array([alpha,beta,gamma]).T
    pos=np.array([x,y,z]).T
    marker=np.array(marker)

    return delta,eta,euler,pos,marker            
        
    

def load_surface(session,tensor_file,sc=2.09,theta_steps=100,phi_steps=50,
                 positive_color=(255,100,100,255),negative_color=(100,100,255,255)):
    
    Delta,Eta,Euler,Pos,Marker=load_tensor(tensor_file)
    
    from chimerax.core.models import Surface
    from chimerax.surface import calculate_vertex_normals,combine_geometry_vntc
    
    geom=list()
    
    for k,(delta,eta,euler,pos,marker) in enumerate(zip(Delta,Eta,Euler,Pos,Marker)):
        if marker==1:
            pc=(100,255,100,255)
            nc=(255,255,100,255)
        else:
            pc=positive_color
            nc=negative_color
        xyz,tri,colors=spherical_surface(\
                                         delta=delta,eta=eta,euler=euler,pos=pos,\
                                         sc=sc,theta_steps=theta_steps,\
                                         phi_steps=phi_steps,\
                                         positive_color=pc,\
                                         negative_color=nc)

        norm_vecs=calculate_vertex_normals(xyz,tri)
        
        geom.append((xyz,norm_vecs,tri,colors))    
        
    xyz,norm_vecs,tri,colors=combine_geometry_vntc(geom)    
    s = Surface('surface',session)
    s.set_geometry(xyz,norm_vecs,tri)
    s.vertex_colors = colors
    session.models.add([s])

    return s


def d2(c=0,s=None,m=None,mp=0):
    """
    Calculates components of the d2 matrix. By default only calculates the components
    starting at m=0 and returns five components, from -2,-1,0,1,2. One may also
    edit the starting component and select a specific final component 
    (mp=None returns all components, whereas mp may be specified between -2 and 2)
    
    d2_m_mp=d2(m,mp,c,s)  #c and s are the cosine and sine of the desired beta angle
    
        or
        
    d2_m_mp=d2(m,mp,beta) #Give the angle directly
    
    Setting mp to None will return all values for mp in a 2D array
    
    (Note that m is the final index)
    """
    
    if s is None:
        c,s=np.cos(c),np.sin(c)
    
    """
    Here we define each of the components as functions. We'll collect these into
    an array, and then call them out with the m and mp indices
    """
    "First, for m=-2"
    
    if m is None or mp is None:
        if m is None and mp is None:
            print('m or mp must be specified')
            return
        elif m is None:
            if mp==-2:
                index=range(0,5)
            elif mp==-1:
                index=range(5,10)
            elif mp==0:
                index=range(10,15)
            elif mp==1:
                index=range(15,20)
            elif mp==2:
                index=range(20,25)
        elif mp is None:
            if m==-2:
                index=range(0,25,5)
            elif m==-1:
                index=range(1,25,5)
            elif m==0:
                index=range(2,25,5)
            elif m==1:
                index=range(3,25,5)
            elif m==2:
                index=range(4,25,5)
    else:
        index=[(mp+2)*5+(m+2)]
    
    out=list()    
    for i in index:
        #mp=-2
        if i==0:x=0.25*(1+c)**2
        if i==1:x=0.5*(1+c)*s
        if i==2:x=np.sqrt(3/8)*s**2
        if i==3:x=0.5*(1-c)*s
        if i==4:x=0.25*(1-c)**2
        #mp=-1
        if i==5:x=-0.5*(1+c)*s
        if i==6:x=c**2-0.5*(1-c)
        if i==7:x=np.sqrt(3/8)*2*c*s
        if i==8:x=0.5*(1+c)-c**2
        if i==9:x=0.5*(1-c)*s
        #mp=0
        if i==10:x=np.sqrt(3/8)*s**2
        if i==11:x=-np.sqrt(3/8)*2*s*c
        if i==12:x=0.5*(3*c**2-1)
        if i==13:x=np.sqrt(3/8)*2*s*c
        if i==14:x=np.sqrt(3/8)*s**2
        #mp=1
        if i==15:x=-0.5*(1-c)*s
        if i==16:x=0.5*(1+c)-c**2
        if i==17:x=-np.sqrt(3/8)*2*s*c
        if i==18:x=c**2-0.5*(1-c)
        if i==19:x=0.5*(1+c)*s
        #mp=2
        if i==20:x=0.25*(1-c)**2
        if i==21:x=-0.5*(1-c)*s
        if i==22:x=np.sqrt(3/8)*s**2
        if i==23:x=-0.5*(1+c)*s
        if i==24:x=0.25*(1+c)**2
        out.append(x)
        
    if m is None or mp is None:
        return np.array(out)
    else:
        return out[0]
