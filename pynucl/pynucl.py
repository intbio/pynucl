# -*- coding: utf-8 -*-
"""
This is the main module of the pynucl package.
Provides entry points for the main class used to analyze nucleosme structures.


Copyright 2020, Alexey Shaytan, Grigoriy Armeev, Moscow Sate Univeristy

"""

#from .base import *
from .nucl_ref_frame_aln import nucl_align
from .nucl_meta_select import create_elem_dict,nucl_sel_expand
from .view_nucl import view_nucl
from .seq_utils import *
from .hist_features import *

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
    def __init__(self, struct, format=None, time=(0,None,1), ref='DNA_align', fullseqs=None, debug=False):
        """
    struct - an MDanalysis universe (single structrue or trajectory)
        can also be: name of a structure file, with 'fromat' - specifying the format if MDanalysis cannot guess.
    ref - reference used for alignment and calculation of deviations
        can be: '1KX5_NRF' - will take the data/1KX5_NRF.pdb - prealigned 1KX5 #TODO - this will be the default in future
        can be: an MD analysis universe - then will use it#TODO
        can be: 'DNA_align' - will use first frame to align it via the nucl_align protocol (3DNA calculations of superhelical axis, etc.)
    fullseqs - information about the full length sequences of the corresponding protein and DNA chains in the structure.
        can be: 'PDB_ID' (e.g. 1KX5) - will try to get info directly from PDB fasta header files (SEQREC records), this might be or might not be the full sequence depending on what the submitters have put into pdb.
        can be: None - just the sequence from the structure will be put instead.
        can be: dict={'segid':'SEQUENCE'}, e.g. {'I':'ATGCATGCATGC...'} or {'A':'GTYHLK...'}
    time - a tuple to specify start, stop and strid along the trajectory to analyze. TODO be careful if to include stop or not(!) should specify here - currently poorly implemented.
        """
        self.time=time
        
        if(isinstance(struct,mda.Universe)):
            self.u=struct
        else:
            self.u=mda.Universe(struct,format=format)
        ###Step 1. Identify components of the system
        # TODO dict={'segid':
        # {'entity':'DNA/RNA/protein/histone/chemical','type':'H4/nuclDNAtop/nuclDNAbottom','variant':'H2A.Z','side':0/1/2}}
        # side 0 - undef, 1 - proximal, 2 - distal.
        # stub here for 1KX5
        self.components={'I':{'entity':'DNA','type':'DNAtop','variant':'alphasat','side':0},
                         'J':{'entity':'DNA','type':'DNAbot','variant':'alphasat','side':0},
                         'A':{'entity':'histone','type':'H3','variant':'canonicalH3','side':1},
                         'E':{'entity':'histone','type':'H3','variant':'canonicalH3','side':2},
                         'B':{'entity':'histone','type':'H4','variant':'canonicalH4','side':1},
                         'F':{'entity':'histone','type':'H4','variant':'canonicalH4','side':2},
                         'C':{'entity':'histone','type':'H2A','variant':'canonicalH2A','side':1},
                         'G':{'entity':'histone','type':'H2A','variant':'canonicalH2A','side':2},
                         'D':{'entity':'histone','type':'H2B','variant':'canonicalH2B','side':1},
                         'H':{'entity':'histone','type':'H2B','variant':'canonicalH2B','side':2}}
    
    
        ###Step 2. Collect sequence information and feature annotation
        self.seqs={}
        self.seqs={k:{} for k in self.components.keys()}#populate with keys
        
        # We should have for every chain:
        # 1) its sequence exactly as in structure 'strseq' - we should extract if from the structure.
        # 2) its sequence corresponding to full protein - currently we will stick to SEQREC from pdb 'fullseq'
        # 3) dictrionary of features for 'fullseq' that can be fed to Pytexshade
        # 4) this dict should also include a feature that shows the 'strseq' on 'fullseq'
        # 5) 'resid_start' - the starting residue number in pdb, i.e. the id of first residue on strseq.
        # 6) 'overhang' - the lenght of left overhang of fullseq over strseq
        #
        # This whole terminalogy is taken from  seqplot/seqplot/pdb_plot.py
        # Where we readjust the data fro plotting as: datafixed['resid']=datafixed['resid']-resid_start+overhang+1-cropseq[0]
        # As a fallback in we cannot get SEQREC - fullseq=strseq
        
        ##Step 2.1. get sequences from the structure:
        # self.seqs={'A':{'seqstr':'AFSDER...'}}
        for k in self.components.keys():
            self.seqs[k]['strseq']=seq_from_mda(self.u,k)
            
        ##Step 2.2. get fullsequences
        if (fullseqs is None): #inherit strseqs
            for k in self.components.keys():
                self.seqs[k]['fullseq']=self.seqs[k]['strseq']
        elif isinstance(fullseqs, str) and len(fullseqs)==4:
            for k in self.components.keys():
                self.seqs[k]['fullseq']=seqrec_from_pdb(fullseqs,k)
        elif isinstance(fullseqs, dict):
            for k in self.components.keys():
                self.seqs[k]['fullseq']=fullseqs[k]
        else:
            logger.error('fullseqs is of unknown type and format')

        ##Step 2.3.
        # Let's calculate overhang of fullseq over strseq and put it to self.seqs[k]['overhang']
        # Let's put the first residue id too.
        for k in self.components.keys():
            self.seqs[k]['overhangL']=get_overhangL(self.seqs[k]['fullseq'],self.seqs[k]['strseq'])
            self.seqs[k]['overhangR']=get_overhangR(self.seqs[k]['fullseq'],self.seqs[k]['strseq'])
            self.seqs[k]['resid_start']=self.u.select_atoms('segid %s and (protein or nucleic)'%k).residues.resids[0]
    
    
        ##Step 2.5. Get sequence features
        # We supply every chain where possible with a list of SeqFeature Objects.
        # These follow the numbering of fullseq
        # 
        self.seq_features={}
        for k in self.components.keys():
            if(self.components[k]['entity']=='histone'):
                self.seq_features[k]=hist_features(self.seqs[k]['fullseq'])
        # TODO: if fullseqs='PDBID', we could grab other features from NCBI
        #    elif(self.components[k]['entity']=='protein'):
        #        self.seq_features[k]=ncbi_features()
        
    
    
    
        ##Step 2.6. Get shading features for use in pytexshade (fullsequence)

        self.shading_features={}
        for k in self.components.keys():
            if(self.components[k]['entity']=='histone'):
                self.shading_features[k]=hist_shade_features(self.seq_features[k])
        # TODO: if fullseqs='PDBID', we could grab other features from NCBI
        #    elif(self.components[k]['entity']=='protein'):
        #        self.seq_features[k]=ncbi_features()
        
        ####Step 3. Make DNA and protein mapping
        #TODO
        #This should create a mapping object it has two dict obj2ref_map ref2obj_map
        #obj2ref_map={(segid,resid):(segid,resid)}
    
        ###Step 4. Make elements dictionary that will be used for meta selections
        #The dict will contain MDanalysis selection strigns for features of every id
        #e.g. alpha1 - will select all alpha1 helices in all histones.
        self.nucl_elements=create_elem_dict(self.components,self.seqs,self.seq_features)
        self.alignment_sel='(%s) and name CA'%self.nucl_elements['hist_folds']
    
        ####Step 5. Place into nucleosome reference frame
        #We have two options here align via calculating DNA superhelix or do it via proxy and align with histone folds. We have in data folder 1KX5_aligned.pdb, but we will need mapping for that.
        #So far will do alignment via DNA #TODO better to optimize it through transformation matrix.
        with tempfile.TemporaryDirectory() as TEMP:
            if(ref=='DNA_align'):
                self.u.select_atoms('all').write(os.path.join(TEMP,'ref.pdb'))
            nucl_align(os.path.join(TEMP,'ref.pdb'),os.path.join(TEMP,'ref_aligned.pdb'),debug=debug)
            refnuc = mda.Universe(os.path.join(TEMP,'ref_aligned.pdb'))
            alignment = align.AlignTraj(self.u, refnuc, select=self.alignment_sel, in_memory=True,start=time[0],stop=time[1],step=time[2])
            alignment.run()      
    
    
    
    def view(self):
        return view_nucl(self.u)
            
class nucltrj(nuclstr):
    """
    Extends nuclstr to import trajectories
    """
    def __init__(self, topol, trj=None, topol_as_first_frame=True, **kwargs):
        if(isinstance(topol,mda.Universe)):
            super().__init__(topol, **kwargs)
        else:
            if(topol_as_first_frame and topol[-3:]=='pdb'):
                opened_trj=mda.Universe(topol,topol,trj)#,in_memory=True
            else:
                opened_trj=mda.Universe(topol,trj)#,in_memory=True
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

