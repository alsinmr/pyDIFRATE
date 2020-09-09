import os
import numpy as np
from chimerax.core.commands import run as rc

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

    delta=np.array(delta)
    eta=np.array(eta)
    euler=np.array([alpha,beta,gamma]).T
    pos=np.array([x,y,z]).T

    return delta,eta,euler,pos            
        
    

def load_surface(session,tensor_file,sc=2.09,theta_steps=100,phi_steps=50,
                 positive_color=(255,100,100,255),negative_color=(100,100,255,255)):
    
    Delta,Eta,Euler,Pos=load_tensor(tensor_file)
    
    from chimerax.core.models import Surface
    from chimerax.surface import calculate_vertex_normals,combine_geometry_vntc
    
    geom=list()
    
    for k,(delta,eta,euler,pos) in enumerate(zip(Delta,Eta,Euler,Pos)):
        xyz,tri,colors=spherical_surface(\
                                         delta=delta,eta=eta,euler=euler,pos=pos,\
                                         sc=sc,theta_steps=theta_steps,\
                                         phi_steps=phi_steps,\
                                         positive_color=positive_color,\
                                         negative_color=negative_color)

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
try:
	load_surface(session,"/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/tensors714808.txt",sc=1,theta_steps=100,phi_steps=50,positive_color=(255, 100, 100, 255),negative_color=(100, 100, 255, 255))
	rc(session,"turn x 90")
	rc(session,"set bgColor white")
	rc(session,"lighting full")
	rc(session,"save /Users/albertsmith/Documents/Dynamics/Lipids/POPC/FinalMD_ana/Ct_orientation_analysis/D2_res83_fr2_t0.1.png ")
except:
	print("Error in chimera script")
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script714808.py")
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/tensors714808.txt")
	rc(session,"exit")
