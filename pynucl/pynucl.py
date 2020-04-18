# -*- coding: utf-8 -*-
"""
This is the main module of the pynucl package.
Provides entry points for the main class used to analyze nucleosme structures.


Copyright 2020, Alexey Shaytan, Grigoriy Armeev, Moscow Sate Univeristy

"""

#from .base import *
from .nucl_ref_frame_aln import nucl_align
from .view_nucl import view_nucl
import MDAnalysis as mda
from MDAnalysis.analysis import align
import pandas as pd
import numpy as np
import tempfile
from io import StringIO
import os
from tqdm import tqdm


import logging
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)s %(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

class nuclstr:
    """
    Main class of the package
    """
    def __init__(self, struct, format=None,time=(None,None),debug=False):
        if(isinstance(struct,mda.Universe)):
            self.u=struct
        else:
            self.u=mda.Universe(struct,format=format)
    
        ####Identify components of the system
        #TODO
    
        ####Make DNA and protein mapping
        #TODO
    
        ####Make elements dictionary
        #TODO
        #stub here
        self.nucl_elements={'hist_folds':'((segid A E and (resid 63:77 or resid 85:114 or resid 120:131)) or (segid B F and (resid 30:41 or resid 49:76 or resid 82:93))  or  (segid C G and (resid 28:39 or resid 59:76 or resid 82:92)) or (segid D H and (resid 34:46 or resid 52:81 or resid 87:99)))'}
        self.alignment_sel='(%s) and name CA'%self.nucl_elements['hist_folds']
    
        ####Place into nucleosome reference frame
        #We have two options here align via calculating DNA superhelix or do it via proxy and align with histone folds. We have in data folder 1KX5_aligned.pdb, but we will need mapping for that.
        #So far will do alignment via DNA #TODO better to optimize it through transformation matrix.
        with tempfile.TemporaryDirectory() as TEMP:
            self.u.select_atoms('all').write(os.path.join(TEMP,'init.pdb'))
            nucl_align(os.path.join(TEMP,'init.pdb'),os.path.join(TEMP,'aligned.pdb'),debug=debug)
            ref = mda.Universe(os.path.join(TEMP,'aligned.pdb'))
            alignment = align.AlignTraj(self.u, ref, select=self.alignment_sel, in_memory=True)
            alignment.run()      
    
    
    
    def view(self):
        return view_nucl(self.u)
            
class nucltrj(nuclstr):
    """
    Extends nuclstr to import trajectories
    """
    def __init__(self, struct, trj, **kwargs):
        if(isinstance(struct,mda.Universe)):
            super().__init__(struct, **kwargs)
        else:
            opened_trj=mda.Universe(struct,trj)
            super().__init__(opened_trj, **kwargs)
#         temp = tempfile.NamedTemporaryFile(suffix='.png')
#         if(debug):
#             print("tempfile created: ",temp.name)

#This is the main enty point of the libary,
#defines classes



# class shadedmsa(object):
#     """Class for storing a shaded image"""

#     def __init__(self, msa,shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False,startnumber=1):
#         temp = tempfile.NamedTemporaryFile(suffix='.png')
#         if(debug):
#             print("tempfile created: ",temp.name)
#         shade.shade_aln2png(msa, filename=temp.name,shading_modes=shading_modes, features=features,title=title,legend=legend, logo=logo,hideseqs=hideseqs,splitN=splitN,setends=setends,ruler=ruler,show_seq_names=show_seq_names,show_seq_length=show_seq_length,funcgroups=funcgroups,rotate=rotate,threshold=threshold,resperline=resperline,margins=margins, density=density,debug=debug,startnumber=startnumber)
#         self.img=open(temp.name, 'rb').read()

#     def _repr_png_(self):
#         return self.img

