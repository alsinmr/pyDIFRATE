#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:21:51 2019

@author: albertsmith
"""

"""
Library of functions to deal with vectors and tensors, used for aligning tensors
and vectors into different frames. We assume all vectors provided are 2D numpy
arrays, with the first dimension being X,Y,Z (we do not deal with time-
dependence in these functions. This is obtained by sweeping over the trajectory.
Frames are processed at each time point separately)
"""


"""
Rotations are 
"""


import numpy as np
from scipy.linalg import svd

#%% Periodic boundary condition check
def pbc_corr(v0,box):
    """
    Corrects for bonds that may be extended across the box. Our assumption is
    that no vector should be longer than half the box. If we find such vectors,
    we will add/subtract the box length in the appropriate dimension(s)
    
    Input should be 3xN vector and 3 element box dimensions
    
    v = pbc_corr(v0,box)
    """
    
    "Copy input, take tranpose for easier calculation"
    v=v0.copy()
    v=v.T
    
    
    i=v>box/2
    ib=np.argwhere(i).T[1]
    v[i]=v[i]-box[ib]
    
    i=v<-box/2
    ib=np.argwhere(i).T[1]
    v[i]=v[i]+box[ib]
        
    return v.T


#%% Periodic boundary condition for positions
def pbc_pos(v0,box):
    """
    Sometimes, we are required to work with an array of positions instead of
    a pair of positions (allowing easy calculation of a vector and determining
    if the vector wraps around the box). In this case, we take differences 
    between positions, and make sure the differences don't yield a step around
    the box edges. The whole molecule, however, may jump around the box after
    this correction. This should matter, since all calculations are orientational,
    so the center position is irrelevant.
    
    Input is a 3xN vector and a 3 element box
    """
    
    v=np.concatenate((np.zeros([3,1]),np.diff(v0,axis=1)),axis=1)
    v=pbc_corr(v,box)
    
    return np.cumsum(v,axis=1)+np.atleast_2d(v0[:,0]).T.repeat(v.shape[1],axis=1)
    
    
#%% Vector normalization
def norm(v0):
    """
    Normalizes a vector to a length of one. Input should be a 3xN vector.
    """
    X,Y,Z=v0
    length=np.sqrt(X**2+Y**2+Z**2)
    
    return v0/length

#%% Reverse rotation direction (passive/active)
def pass2act(cA,sA,cB,sB=None,cG=None,sG=None):
    """
    After determining a set of euler angles, we often want to apply them to go
    into the reference frame corresponding to those angles. This requires
    reversing the rotation, performed by this function
    
    -gamma,-beta,-alph=pass2act(alpha,beta,gamma)
    
    or 
    
    cG,-sG,cG,-sB,cA,-sA=pass2act(cA,sA,cB,sB,cG,sG)
    """
    
    if sB is None:
        return -cB,-sA,-cA
    else:
        return cG,-sG,cB,-sB,cA,-sA
    
#%% Frame calculations
def getFrame(v1,v2=None,return_angles=False):
    """
    Calculates the sines and cosines of the euler angles for the principle axis
    system of a frame defined by one or two vectors. The new frame has v1 along
    its z-axis and if a second vector is provided, then the second vector lies
    in the xz plane of the frame.
    
    We use zyz convention (alpha,beta,gamma), where rotation into a frame is 
    achieved by first applying gamma:
        X,Y=cos(gamma)*X+sin(gamma)*Y,-sin(gamma)*X+cos(gamma)*Y
    Then applying beta:
        X,Z=cos(beta)*X-sin(beta)*Z,sin(beta)*X+cos(beta)*Z
    Finally alpha:
        X,Y=cos(alpha)*X+sin(alpha)*Y,-sin(alpha)*X+cos(alpha)*Y
        
    gamma=arctan2(Y,X)
    beta=arccos(Z)
    alpha=arctan(Y1,X1) (Y1 and X1 after applying gamma and beta!)
    
    Note that we do not return alpha,beta,gamma! Instead, we return 
    cos(alpha),sin(alpha),cos(beta),sin(beta),cos(gamma),sin(gamma)!
    
    If only one vector is provided, then we simply require that this 
    vector lies along z, achieved by rotating the shortest distance to the 
    z-axis. Then, the euler angles are (-gamma,beta,gamma)
    
    
    cA,sA,cB,sB,cG,sG = getFrame(v1,v2)
    
        or
        
    cG,-sG,cB,sB,cG,sG = getFrame(v1)
    
    Finally, if you need the angles themselves:
        
    alpha,beta,gamma = getFrame(v1,v2,return_angles=True)
    """
    
    "Set nan values along z (these indicate an unused frame)"
    ii=np.isnan(v1[0,:,0])
        
    "Normalize"
    X,Y,Z=norm(v1)
    
    "Gamma"
    lenXY=np.sqrt(X**2+Y**2)
    i=lenXY==0
    lenXY[i]=1  #cG and sG will be 0
    cG,sG=-X/lenXY,Y/lenXY
    cG[i]=1. #Set cG to 1 where cG/sG is undefined (gamma=0)
    
    "Beta"
    cB,sB=Z,np.sqrt(1-Z**2)
    
    "Alpha"
    if v2 is None:
        cA,sA=cG,-sG #If only one vector, alpha=-gamma
    else:
        v2=Rz(v2,cG,-sG)
        X,Y,_=Ry(v2,cB,-sB)
        
        lenXY=np.sqrt(X**2+Y**2)
        i=lenXY==0
        lenXY[i]=1  #cA and sA will still be 0
        cA,sA=X/lenXY,-Y/lenXY
        cA[i]=1.
        
#    cA[ii]=1
#    sA[ii]=0
#    cB[ii]=1
#    sB[ii]=0
#    cG[ii]=1
#    sG[ii]=0
    if return_angles:
        return np.arctan2(sA,cA),np.arctan2(sB,cB),np.arctan2(sG,cG)
    else:
        return cA,sA,cB,sB,cG,sG
    

#%% Apply/invert rotations     
def Rz(v0,c,s=None):
    """
    Rotates a vector around the z-axis. One must provide the vector(s) and either
    the angle itself, or the cosine(s) and sine(s) of the angle(s). The number
    of vectors must match the number of angles, or only one angle is provided
    
    v=Rz(v0,c,s)
    
        or
        
    v=Rz(v0,theta)
    """
    
    if s is None:
        c,s=np.cos(c),np.sin(c)
        
    X,Y,Z=v0.copy()
    
    X,Y=c*X+s*Y,-s*X+c*Y
    
    return np.array([X,Y,Z])

def Ry(v0,c,s=None):
    """
    Rotates a vector around the y-axis. One must provide the vector(s) and either
    the angle itself, or the cosine(s) and sine(s) of the angle(s). The number
    of vectors must match the number of angles, or only one angle is provided
    
    v=Ry(v0,c,s)
    
        or
        
    v=Ry(v0,theta)
    """
    
    if s is None:
        c,s=np.cos(c),np.sin(c)
        
    X,Y,Z=v0.copy()
    
    X,Z=c*X-s*Z,s*X+c*Z
    
    return np.array([X,Y,Z])

def R(v0,cA,sA,cB,sB=None,cG=None,sG=None):
    """
    Rotates a vector using ZYZ convention. One must provide the vector(s) and 
    either the euler angles, or the cosine(s) and sine(s) of the angle(s). The 
    number of vectors must match the number of angles, or only one angle is 
    provided for alpha,beta,gamma (or the sines/cosines of alpha,beta,gamma)
    
    v=R(v0,cA,sA,cB,sB,cG,sG)
    
        or
        
    v=R(v0,alpha,beta,gamma)
    """
    if v0 is None:
        return None
    
    if sB is None:
        cA,sA,cB,sB,cG,sG=np.cos(cA),np.sin(cA),np.cos(sA),np.sin(sA),np.cos(cB),np.sin(cB)
        
    v=Rz(v0,cA,sA)
    v=Ry(v,cB,sB)
    v=Rz(v,cG,sG)
    
    return v


def R2euler(R,return_angles=False):
    """
    Input a rotation matrix in cartesian coordinates, and return either the
    euler angles themselves or their cosines and sines(default)
    
    cA,sA,cB,sB,cG,sG = R2euler(R)
    
        or
    
    alpha,beta,gamma = R2euler(R,return_angles=True)
    
    R can be a list of matrices
    """
    
    R = np.array([R]) if np.ndim(R)==2 else np.array(R)
    
    
    """
    Note that R may be the result of an eigenvector decomposition, and does
    not guarantee that R is a proper rotation matrix. We can check the sign
    on the determinant: if it is 1, it's a proper rotation, if it's -1, it's not
    Then, we just multiply each matrix by the result to have only proper
    rotations.
    
    Currently, this correction is being made elsewhere. Uncomment if it should
    be applied here
    """
    sgn=np.sign(np.linalg.det(R))
    
    R=(R.T*sgn).T
    
    cB=R[:,2,2]
    cB[cB>1]=1.
    cB[cB<-1]=-1.
    sB=np.sqrt(1.-cB**2)
    i=sB!=0
    cA,sA,cG,sG=[np.ones(i.shape),np.zeros(i.shape),np.ones(i.shape),np.zeros(i.shape)]
    cA[i]=R[i,2,0]/sB[i]
    sA[i]=R[i,2,1]/sB[i]
    cG[i]=-R[i,0,2]/sB[i]
    sG[i]=R[i,1,2]/sB[i]
    cG[np.logical_not(i)]=R[np.logical_not(i),0,0]
    sG[np.logical_not(i)]=R[np.logical_not(i),0,1]
    
    if return_angles:
        return np.array((np.arctan2(sA,cA),np.arctan2(sB,cB),np.arctan2(sG,cG)))
    else:
        return np.array((cA,sA,cB,sB,cG,sG))
    
def R2vec(R):
    """
    Given a rotation matrix, R, this function returns two vectors, v1, and v2
    that have been rotated from v10=[0,0,1] and v20=[0,1,0]
    
    v1=np.dot(R,v10)
    v2=np.dot(R,v20)
    
    If a frame is defined by a rotation matrix, instead of directly by a set of
    vectors, then v1 and v2 have the same Euler angles to rotate back to their
    PAS as the rotation matrix
    
    R may be a list of rotation matrices
    
    Note: v1, v2 are trivially given by R[:,:,2] and R[:,:,0]
    """
    R = np.array([R]) if np.ndim(R)==2 else np.array(R)
    
    v1=R[:,:,2]
    v2=R[:,:,0]
    
    return v1.T,v2.T
    
    
#%% Tensor operations
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
    if mp==-2:
        def d2m2(c,s):return 0.25*(1+c)**2
        def d2m1(c,s):return 0.5*(1+c)*s
        def d20(c,s):return np.sqrt(3/8)*s**2
        def d21(c,s):return 0.5*(1-c)*s
        def d22(c,s):return 0.25*(1-c)**2
    elif mp==-1:
        def d2m2(c,s):return -0.5*(1+c)*s
        def d2m1(c,s):return c**2-0.5*(1-c)
        def d20(c,s):return np.sqrt(3/8)*2*c*s
        def d21(c,s):return 0.5*(1+c)-c**2
        def d22(c,s):return 0.5*(1-c)*s
    elif mp==0:
        def d2m2(c,s):return np.sqrt(3/8)*s**2
        def d2m1(c,s):return -np.sqrt(3/8)*2*s*c
        def d20(c,s):return 0.5*(3*c**2-1)
        def d21(c,s):return np.sqrt(3/8)*2*s*c
        def d22(c,s):return np.sqrt(3/8)*s**2
    elif mp==1:
        def d2m2(c,s):return -0.5*(1-c)*s
        def d2m1(c,s):return 0.5*(1+c)-c**2
        def d20(c,s):return -np.sqrt(3/8)*2*s*c
        def d21(c,s):return c**2-0.5*(1-c)
        def d22(c,s):return 0.5*(1+c)*s
    elif mp==2:
        def d2m2(c,s):return 0.25*(1-c)**2
        def d2m1(c,s):return -0.5*(1-c)*s
        def d20(c,s):return np.sqrt(3/8)*s**2
        def d21(c,s):return -0.5*(1+c)*s
        def d22(c,s):return 0.25*(1+c)**2
        
    if m is None:
        return np.array([d2m2(c,s),d2m1(c,s),d20(c,s),d21(c,s),d22(c,s)])
    elif m==-2:
        return d2m2(c,s)
    elif m==-1:
        return d2m1(c,s)
    elif m==0:
        return d20(c,s)
    elif m==1:
        return d21(c,s)
    elif m==2:
        return d22(c,s)
    
def D2(cA=0,sA=0,cB=0,sB=None,cG=None,sG=None,m=None,mp=0):
    """
    Calculates components of the Wigner rotation matrix from Euler angles or
    from the list of sines and cosines of those euler angles. All vectors must
    be the same size (or have only a single element)
    
    mp and m should be specified. m may be set to None, so that all components
    are returned in a 2D array
    
    D2_m_mp=D2(m,mp,cA,sA,cB,sB,cG,sG)  #Provide sines and cosines
    
        or
        
    D2_m_mp=D2(m,mp,alpha,beta,gamma) #Give the angles directly
    
    (Note that m is the final index)
    """
    
    if sB is None:
        cA,sA,cB,sB,cG,sG=np.cos(cA),np.sin(cA),np.cos(sA),np.sin(sA),np.cos(cB),np.sin(cB)
        
    d2c=d2(cB,sB,m,mp)
    
    "Rotation around z with alpha (mp)"
    if mp!=0:
        ea=cA-1j*np.sign(mp)*sA
        if np.abs(mp)==2:
            ea=ea**2
    else:
        ea=1

    "Rotation around z with gamma (m)"
    if m is None:
        eg1=cG-1j*sG
        egm1=cG+1j*sG
        eg2=eg1**2
        egm2=egm1**2
        eg0=np.ones(eg1.shape)
        eg=np.array([egm2,egm1,eg0,eg1,eg2])
     
    return ea*d2c*eg
    

def D2vec(v1,v2=None,m=None,mp=0):
    """
    Calculates the Wigner rotation elements that bring a vector or vectors from
    their own principle axis system into a reference frame (whichever frame
    v1 and v2 are defined in)
    """
    
    cA,sA,cB,sB,cG,sG=getFrame(v1,v2)
    "I think these are already the passive angles above"
    
    return D2(cA,sA,cB,sB,cG,sG,m,mp)

def Spher2Cart(rho):
    """
    Takes a set of components of a spherical tensor and calculates its cartesian
    representation (as a vector, with components in order of Axx,Axy,Axz,Ayy,Ayz)
    
    Input may be a list (or 2D array), with each new column a new tensor
    """
    
    rho=np.array(rho,dtype=complex)

    M=np.array([[0.5,0,-np.sqrt(1/6),0,0.5],
                 [0.5*1j,0,0,0,-0.5*1j],
                 [0,0.5,0,-0.5,0],
                 [-0.5,0,-np.sqrt(1/6),0,-.5],
                 [0,.5*1j,0,.5*1j,0]])
    return np.dot(M,rho).real
    
    
def Spher2pars(rho,return_angles=False):
    """
    Takes a set of components of a spherical tensor and calculates the parameters
    describing that tensor (delta,eta,alpha,beta,gamma)
    
    
    delta,eta,cA,sA,cB,sB,cG,sG=Spher2pars(rho)
    
        or
        
    delta,eta,alpha,beta,gamma=Spher2pars(rho,return_angles=True)
    
    
    Input may be a list (or 2D array), with each new column a new tensor
    
    """

    A0=np.atleast_2d(Spher2Cart(rho)) #Get the cartesian tensor

    R=list()
    delta=list()
    eta=list()
    
    
    for k,x in enumerate(A0.T):
        Axx,Axy,Axz,Ayy,Ayz=x
        A=np.array([[Axx,Axy,Axz],[Axy,Ayy,Ayz],[Axz,Ayz,-Axx-Ayy]])    #Full matrix
        D,V=np.linalg.eigh(A)   #Get eigenvalues, eigenvectors 
        i=np.argsort(np.abs(D))
        D,V=D[i[[1,0,2]]],V[:,i[[1,0,2]]]     #Ordering is |azz|>=|axx|>=|ayy|
        "V should have a determinant of +1 (proper vs. improper rotation)"
        V=V*np.sign(np.linalg.det(V))
        delta.append(D[2])
        eta.append((D[1]-D[0])/D[2])
        R.append(V)
    
    delta=np.array(delta)
    eta=np.array(eta)
    euler=R2euler(R)
    
    if return_angles:
        cA,sA,cB,sB,cG,sG=euler
        euler=np.array([np.arctan2(sA,cA),np.arctan2(sB,cB),np.arctan2(sG,cG)])
       
    return np.concatenate(([delta],[eta],euler),axis=0)
        

#%% RMS alignment
def RMSalign(v0,vref):
    """
    Returns the optimal rotation matrix to rotate a set of vectors v0 to a set 
    of reference vectors, vref
    
    R=alignR(v0,vref)
    
    Uses the Kabsch algorithm. Assumes *vectors*, with origins at zero, not 
    points, so that no translation will be performed
    (reference https://en.wikipedia.org/wiki/Kabsch_algorithm)
    
    We minimize
    np.sum((np.dot(R,v0.T).T-vref)**2)
    """
    
    H=np.dot(v0,vref.T)
    
    U,S,Vt=svd(H)
    V=Vt.T
    Ut=U.T
    
    d=np.linalg.det(np.dot(V,Ut))
    m=np.eye(3)
    m[2,2]=d    #This is supposed to ensure a right-handed coordinate system
                #But I think it could equivalently be thought of as making this a proper rotation(??)
    
    R=np.dot(V,np.dot(m,Ut))
    return R


#%% Fit points to a plane
def RMSplane(v,weight=None):
    """
    For a set of points (v: 3xN array), calculates the normal vector for a plane
    fitted to that set of points. May include a weighting (weight: N elements)
    """
    v=np.array(norm(v))
    
    "Default, uniform weighting"
    if weight is None:
        weight=np.ones(v.shape[1])
        
    "Subtract away the centroid"
    v=(v.T-v.mean(axis=1)).T
    
    """Applying weighting, taking singular value decomposition, return
    row of U corresponding to the smallest(last) singular value"""
    return svd(v*weight)[0].T[2]

#%% Get principle axes of moment of inertia
def principle_axis_MOI(v):
    """
    Calculates the principle axis system of the moment of inertia, without
    considering weights of individual particles. A 3xN numpy array should be 
    provided. The smallest component of the moment of inertia is returned in the
    0 element, and largest in the 2 element
    
    Note- the directions of the principle axes can switch directions (180 deg)
    between frames, due to the symmetry of the MOI tensor. This can be corrected
    for with a reference vector. The dot product of the reference vector and the
    vector for a given frame should remain positive. If it doesn't, then switch
    the direction of the vector (v=v*np.sign(np.dot(v.T,v)))
    """
    
    """
    Ixx=sum_i m_i*(y_i^2+z_i^2)
    Iyy=sum_i m_i*(x_i^2+z_i^2)
    Izz=sum_i m_i*(x_i^2+y_i^2)
    Ixy=Iyx=-sum_i m_i*x_i*y_i
    Ixz=Izx=-sum_i m_i*x_i*z_i
    Iyz=Izy=-sum_i m_i*y_i*z_i
    """
    
    
    v=v-np.atleast_2d(v.mean(axis=1)).T.repeat(v.shape[1],axis=1) #v after subtracting center of mass
    
    H=np.dot(v,v.T)
    
    I=-1*H
    I[0,0]=H[1,1]+H[2,2]
    I[1,1]=H[0,0]+H[2,2]
    I[2,2]=H[1,1]+H[0,0]
    _,V=np.linalg.eigh(I)
    
    return V

#%% Project onto axis
def projZ(v0,vr=[0,0,1]):
    """
    Takes the projection of a vector, v0, onto another vector, vr.
    
    Input should be 3xN vectors (vnorm can also be a 1D, 3 element vector, or
    both inputs can be 3 element vectors).
    
    Input does not need to be normalized, but also note that output is not 
    normalized
    
    Default project is along z
    """
    v0=np.atleast_2d(v0)
    if np.ndim(vr)==1 or np.shape(vr)[1]==1:    #Make matrices the same size
        vr=np.atleast_2d(vr).T.repeat(np.shape(v0)[1],axis=1) 
    
    return np.atleast_2d((norm(v0)*norm(vr)).sum(axis=0))*vr

#%% Project onto plane
def projXY(v0,vnorm=[0,0,1]):
    """
    Takes the projection of a vector, v0, onto a plane defined by its normal 
    vector, vnorm. 
    
    Input should be 3xN vectors (vnorm can also be a 1D, 3 element vector, or
    both inputs can be 3 element vectors).
    
    Input does not need to be normalized, but also note that output is not 
    normalized
    """
    
    return v0-projZ(v0,vnorm)  
    
#%% Sort by distance
def sort_by_dist(v,maxd=1e4):
    """
    Returns an index that sorts a set of points such that each point is next
    to its nearest neighbors in space in the vector itself (although points will
    not repeat, so this has exceptions)
    
    Searchest for the point closest to (-Inf,-Inf,-Inf), then looks for its nearest
    neighbor, and then searchest for the nearest neighhbor of the next point
    (etc...)
    
    The purpose is that we can take a set of points, and take the difference in
    position of each one to generate a set of vectors (which may be subsequently
    corrected for periodic boundary conditions)
    
    Returns the sorting index, as oppposed to the vector itself
    
    i=sort_by_dist(v)
    
    such that v_sorted=v[i]
    
    Note 1: not a highly efficient algorithm- recommended only for setup of a 
    vector calculation, but should avoided inside loop over a trajectory

    Note 2: We presume here that the dimensions can't be larger than 1e4, rather
    than assuming the dimension is arbitrarily large (creating some numeric
    problems). If this is not true, set maxd to an appropriate value
    """
    
    v=v.copy()  #Avoid editing the vector itself...never quite sure when this is necessary
    X,Y,Z=v.T
    
    ref=maxd*2
    
    i=list()
    "Find the most negative element"
    i.append(np.argmin((X+ref)**2+(Y+ref)**2+(Z+ref)**2))
    
    for _ in range(X.size-1):
        x,y,z=X[i[-1]].copy(),Y[i[-1]].copy(),Z[i[-1]].copy()
        "Set the currect vector to be far away"
        X[i[-1]],Y[i[-1]],Z[i[-1]]=ref*np.array([-2,-2,-2])
        "Find the nearest vector"
        i.append(np.argmin((X-x)**2+(Y-y)**2+(Z-z)**2))
        
    return i
    
    
    
    
    