#Functions to align nucleosome to reference frame by finding a nucleosome axis
#via 3DNA and MDAnalysis

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np
import os
import re
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.optimize import minimize
from .transformations import euler_matrix
from .transformations import translation_matrix
from .transformations import rotation_matrix
from numpy import random
import os
import logging
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)s %(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

from importlib import reload 
import DNAtools

class Noconv(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, message):
        self.message = message
        
def file_len(fname):
    for i, l in enumerate(fname):
        pass
    fname.seek(0)
    return i + 1

def nucl_align(sys_ref_pdb,sys_ref_pdb_aligned,fp_ref=None,x3dna_fp_out=None,init_x0=[0.5]*5,attempts=10,debug=False):
    x0=init_x0
    for i in range(attempts):
        d=nucl_align_worker(sys_ref_pdb,sys_ref_pdb_aligned,fp_ref=fp_ref,x3dna_fp_out=x3dna_fp_out,init_x0=x0,debug=debug)
        if(d<1000):
            ref_fr=mda.Universe(sys_ref_pdb_aligned,in_memory=True)
            sel=ref_fr.select_atoms("nucleic")
            if(sel.positions[0][2])<0:
                if i==(attempts-1):
                    raise Noconv("Error: could not converge to position with DNA end z>0 in %d attemps!!!"%attempts)
                else:
                    logger.debug("\n==DNA start is at z<0, we want it to be at z>0, trying to reinitialize with random numbers==\n")
                    x0=random.rand(5)*5.0
                    logger.debug(x0)
                    logger.debug("Retrying...")
            else:
                break
        else:
            if i==(attempts-1):
                raise Noconv("Error: did not converge after %d attemps!!!"%attempts)
            else:
                logger.debug("\n==Did not converge, trying to reinitialize with random numbers==\n")
                x0=random.rand(5)*5.0
                logger.debug(x0)
                logger.debug("Retrying...")
        
            


def nucl_align_worker(sys_ref_pdb,sys_ref_pdb_aligned,fp_ref=None,x3dna_fp_out=None,init_x0=[0.5]*5,debug=False):
    """
    sys_ref_pdb - PDB file of nucleosome
    sys_ref_pdb_aligned - file to save aligned PDB file of nucleosome
    fp_ref - pdb to get reference base pairing from.
    """
    # find center of DNA molecule
    ref_fr=mda.Universe(sys_ref_pdb,in_memory=True)
    sel=ref_fr.select_atoms("nucleic")
    logger.debug("DNA chains consist of %d "%sel.n_residues + "residues")
    logger.debug("DNA consists of %d"%(sel.n_residues/2) + "basepairs")
    cent_bp_index=int(((sel.n_residues/2)-1)/2) #zero-based, if 146, will get 72, which is ok for 1aoi, but should check for others.
    logger.debug("Cetner bp index %d"%cent_bp_index)
    if(fp_ref):
        logger.debug("Using %s to extract reference base pairing"%fp_ref)
        ref=DNAtools.X3DNA_find_pair(fp_ref)
    else:
        if(x3dna_fp_out):
            ref=x3dna_fp_out
        else:
            logger.debug("Calculating reference base pairing from original file")
            ref=DNAtools.X3DNA_find_pair(sys_ref_pdb)
   
    logger.debug('\n===BP identification file is %s\n'%ref)
    logger.debug(ref.read())
    ref.seek(0)
    #Check if the number of identified base pairs is equal to the real number of base pairs
    ident_bp=file_len(ref)-8
    logger.debug("identified bp = %d"%ident_bp)
    if not (ident_bp==sel.n_residues/2):
        raise AssertionError("The number of base pairs identified by X3DNA is not equal to that expected from DNA length:%d,%d"%(ident_bp,sel.n_residues/2))
    
    df=DNAtools.X3DNA_analyze(sys_ref_pdb,ref)
    #We need to check that 3DNA identified all the base pairs.
    #NOTE: by default this is not the case for 1kx5 structure. we tweak the 3DNA parameters in it config files.
    assert(len(df)==sel.n_residues/2)
    # Generate transformation matrix 


    def calc_dist(v):
        """
        Calculates sum(d-<d>)^2 from Z-axis
        """
        np=v.shape[1]
        #Get average
        dav=0
        for i in range(np):
            dav=dav+(v[0,i]**2+v[1,i]**2)**0.5
        dav=dav/np
        dsum=0
        for i in range(np):
            dsum=dsum+((v[0,i]**2+v[1,i]**2)**0.5-dav)**2
        return dsum



    #we need rotation and calculation
    def fun(p):
        """
        does rotation and calls calc_dist
        p[0] - alpha, p[1] - beta, p[2] -gamma, p[3] - x transl, p[4] - y transl
        """
        Re = euler_matrix(p[0], p[1], p[2], 'rzxz')
        T=translation_matrix([p[3],p[4],0])
        M=np.dot(T,Re)
        vn=np.dot(M,np.concatenate((vertices,np.array([[1]*vertices.shape[1]])),0))[:-1,:]
        return(calc_dist(vn))

    def funv(p):
        """
        does rotation 
        p[0] - alpha, p[1] - beta, p[2] -gamma, p[3] - x transl, p[4] - y transl
        """
        Re = euler_matrix(p[0], p[1], p[2], 'rzxz')
        T=translation_matrix([p[3],p[4],0])
        M=np.dot(T,Re)
        vn=np.dot(M,np.concatenate((vertices,np.array([[1]*vertices.shape[1]])),0))[:-1,:]
        return(vn)

    def funm(p):
        """
        returns rotation matrix
        p[0] - alpha, p[1] - beta, p[2] -gamma, p[3] - x transl, p[4] - y transl
        """
        Re = euler_matrix(p[0], p[1], p[2], 'rzxz')
        T=translation_matrix([p[3],p[4],0])
        M=np.dot(T,Re)
        return(M)


    vertices=(df[['x','y','z']].values).T
    logger.debug("Initial bp positions: "+vertices.__repr__())
    #Now we need to get nucleosome axis
    logger.debug("Intial distance var: %f"%calc_dist(vertices))
    x0 = init_x0
    
    res = minimize(fun, x0,method='Nelder-Mead',options={'maxiter':10000})
#     res = minimize(fun, x0,method='Nelder-Mead')
#     res = minimize(fun, x0,method='Powell')

#     print(res)
#     print(res.x)
    logger.debug("Distance var after minim: %f"%fun(res.x))


 #   if(fun(res.x)>1000):
 #       raise Noconv("The value of dist var more than 1000 generally means we did not find the minimum, try changing x0!!!")


    vnew=funv(res.x)
    MAT1=funm(res.x)
    #print("Bp positions after first transformation:",vnew)



    #Now let's get a perpendicular

    cbp=vnew[:,int(cent_bp_index)]
    z=cbp[2]
    y=cbp[1]
    x=cbp[0]
    #print("Center bp coord:",x,y,z)
    #we need this point to lie on Y axis! According to Zhurkin

    #Seems that it is better to take
    #Z as geometrical center of all BP
    z=np.mean(vnew[2,:])
    #print("Center Z coord:",z)
    tmz=translation_matrix([0,0,-z])
    MAT2=np.dot(tmz,MAT1)
    logger.debug('%f'%x)
    logger.debug('%f'%y)
    rmzt=rotation_matrix(3.14159/2-np.angle(x+y*1j),[0,0,1])
    # rmzt=rotation_matrix(3.14,[0,0,1])

    MAT3=np.dot(rmzt,MAT2)
    #print("Final transformation matrix\n",MAT3)
    # MAT3=MAT2
    # t=(-0.9996969699859619, -0.002169886138290167, 0.024519488215446472, 0.0, -0.0014522717101499438, 0.9995711445808411, 0.0292470995336771, 0.0, -0.02457243576645851, 0.029202627018094063, -0.9992714524269104, 0.0, 0.3103618621826172, -0.13102436065673828, -0.9151911735534668, 1.0)
    #MATFIT=(MAT3.T).reshape((1,16))[0] #Magic
    #print(MATFIT) # this is the final matix
    #Let's try to rotate using MDanalysis
    sel=ref_fr.select_atoms("all")
    num_atoms=sel.positions.shape[0]
    d4_vec=np.append(sel.positions.T,np.ones((1,num_atoms)),axis=0)
    logger.debug(d4_vec.__repr__())
    d4_vec_new=np.dot(MAT3,d4_vec)
    logger.debug(d4_vec_new.__repr__())
    sel.positions=d4_vec_new[0:3].T
    # Set final PDB-file name 
    sel.write(sys_ref_pdb_aligned)
    return(fun(res.x))



