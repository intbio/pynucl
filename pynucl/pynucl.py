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
from .entity_detector import detect_entities,check_nucleosomes_sides,get_base_pairs_dict
from .utils.structure_constants import IonNames

import MDAnalysis as mda
from MDAnalysis.analysis import align
import pandas as pd
import numpy as np
import tempfile
from io import StringIO
import os
from tqdm.auto import tqdm
import uuid
import socket
import pkg_resources




import logging
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)s %(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

class nuclstr:
    """
    Main class of the package
    """
    def __init__(self, struct,name=None, format=None, time=(0,None,1), ref='1KX5_NRF', fullseqs=None, debug=False,skipaln=False,aln_sel=None,auto_detect_entities=False,check_dyad=True):
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
        if not hasattr(self, 'name'):
            self.name=name
        
        self.time=time
        
        if(isinstance(struct,mda.Universe)):
            self.u=struct
        else:
            self.u=mda.Universe(struct,format=format)
        self.frames=len(self.u.trajectory)
        print('%d frames loaded for %s'%(self.frames, self.name if hasattr(self,'name') else 'system'))
        ###Step 1. Identify components of the system
        # TODO dict={'segid':
        # {'entity':'DNA/RNA/protein/histone/chemical','type':'H4/nuclDNAtop/nuclDNAbottom','variant':'H2A.Z','side':0/1/2}}
        # side 0 - undef, 1 - proximal, 2 - distal.
        # stub here for 1KX5
        if auto_detect_entities:
            self.components=detect_entities(self.u)
        else:
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
        selection=self.u.select_atoms(f'(protein or nucleic) and (not resname {" ".join(IonNames)})')
        for k in self.components.keys():
            if(self.components[k]['entity'] in ['histone','DNA']):
                segment=selection.select_atoms(f'segid {k}')
                segment=segment[(segment.altLocs=='') | (segment.altLocs=='A')]
                strseq=seq_from_mda(segment,k)
                resids=segment.residues.resids
                self.seqs[k]['strseq']=strseq[0]+''.join(['X'*(diff-1) + strseq[i+1] for i,diff in enumerate(np.diff(resids))])

        ##Step 2.2. get fullsequences
        if (fullseqs is None): #inherit strseqs
            for k in self.components.keys():
                if(self.components[k]['entity'] in ['histone','DNA']):
                    #self.seqs[k]['fullseq']=self.seqs[k]['strseq']
                    self.seqs[k]['fullseq']=self.seqs[k]['strseq']
                    self.seqs[k]['overhangL']=self.seqs[k]['overhangR']=0
                    self.seqs[k]['resid_start']=self.u.select_atoms('segid %s and (protein or nucleic)'%k).residues.resids[0]
        else:
            if isinstance(fullseqs, str) and len(fullseqs)==4:
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
                if self.components[k]['type']!='H1':
                    self.seq_features[k]=hist_features(self.seqs[k]['fullseq'])
        # TODO: if fullseqs='PDBID', we could grab other features from NCBI
        #    elif(self.components[k]['entity']=='protein'):
        #        self.seq_features[k]=ncbi_features()        
        
        # last check: need to detect sides and NLPs (nucleosome like particles)
        if detect_entities:
            # this element dict is not valid later as all histones in it are on side 1
            elements_dict=create_elem_dict(self.components,self.seqs,self.seq_features)
            self.components=check_nucleosomes_sides(self.components,elements_dict,self.u,self.seq_features)
    
    
    
        ##Step 2.6. Get shading features for use in pytexshade (fullsequence)

        self.shading_features={}
        for k in self.components.keys():
            if(self.components[k]['entity']=='histone'):
                if self.components[k]['type']!='H1':
                    self.shading_features[k]=hist_shade_features(self.seq_features[k])
                else:
                    self.shading_features[k]=[]
            else:
                self.shading_features[k]=[]
        # TODO: if fullseqs='PDBID', we could grab other features from NCBI
        #    elif(self.components[k]['entity']=='protein'):
        #        self.seq_features[k]=ncbi_features()
        
        ####Step 3. Make DNA and protein mapping
        #TODO - do the same for DNA in sequence space if possible(!)
        #Let's start with protein mapping of histones.
        #This should create a mapping object it has two dict map_obj2ref_resids={(segid,resid):(segid,resid)}
        #map_ref2obj_resids={(segid,resid):(segid,resid)}
        #map_obj2ref_fullseq={(segid,0-based pos):(segid,0-based pos)} , map_obj2ref_fullseq
        
        #Let's start with mapping fullseq to 1kx5 seqs (the latter are identical to reall fullseqs of Xenopus except for H2B - we intentionally retain the sequences that are in 1kx5 so that we could easier map them to resids later)
        map_1kx5={('H3',1):'A',('H4',1):'B',('H2A',1):'C',('H2B',1):'D',('H3',2):'E',('H4',2):'F',('H2A',2):'G',('H2B',2):'H'}
        
        self.map_obj2ref_fullseq={}
        self.map_ref2obj_fullseq={}
        self.map_obj2ref_resids={}
        self.map_ref2obj_resids={}
        for k,v in self.components.items():
            if v['entity']=='histone':
                if v['type']!='H1':
                    self.map_obj2ref_fullseq[k]= map_seqs(self.seqs[k]['fullseq'],ss_templ_dict[v['type']]['seq'],k,map_1kx5[(v['type'],v['side'])])
                    self.map_obj2ref_resids[k]={(k2[0], k2[1]+self.seqs[k]['resid_start']-self.seqs[k]['overhangL']): (v2[0],v2[1]+1) for k2, v2 in self.map_obj2ref_fullseq[k].items()}
                    self.map_ref2obj_fullseq[k] = {v2: k2 for k2, v2 in self.map_obj2ref_fullseq[k].items()}
                    self.map_ref2obj_resids[k] = {v2: k2 for k2, v2 in self.map_obj2ref_resids[k].items()}
        
    
        ###Step 4. Make elements dictionary that will be used for meta selections
        #The dict will contain MDanalysis selection strigns for features of every id
        #e.g. alpha1 - will select all alpha1 helices in all histones.
        
        #These are adjusted internally to switch from
        self.nucl_elements=create_elem_dict(self.components,self.seqs,self.seq_features)
        
        if(aln_sel is None):
            self.alignment_sel=nucl_sel_expand('(alpha1 or alpha2 or alpha3) and name CA',self.nucl_elements)
        else:
            self.alignment_sel=nucl_sel_expand(aln_sel,self.nucl_elements)
    
        ####Step 5. Place into nucleosome reference frame
        #We have two options here align via calculating DNA superhelix or do it via proxy and align with histone folds. We have in data folder 1KX5_aligned.pdb, but we will need mapping for that.
        #So far will do alignment via DNA #TODO better to optimize it through transformation matrix.
        if(not skipaln):
            with tempfile.TemporaryDirectory() as TEMP:
                if(ref=='DNA_align'):
                    self.u.select_atoms('all').write(os.path.join(TEMP,'ref.pdb'))
                    nucl_align(os.path.join(TEMP,'ref.pdb'),os.path.join(TEMP,'ref_aligned.pdb'),debug=debug)
                    self.alignment_sel_dict=self.alignment_sel
                    refnuc = mda.Universe(os.path.join(TEMP,'ref_aligned.pdb'))
    #             os.system('cp '+TEMP+'/ref_aligned.pdb'+' ref_aligned.pdb')
                if(ref=='1KX5_NRF'): # we will align to 1KX5_NRF available in data folder
#                     pass
                    DATA_PATH = pkg_resources.resource_filename('pynucl', 'data/')
                    refnuc = mda.Universe(os.path.join(DATA_PATH,'1KX5_NRF.pdb'))
#                     align.fasta2select('sequences.aln')
                    #we will need to provide for aligment a selection in 1KX5 and analogous selection in the provided structure
                    #we take self.alignment_sel - from that sel we need to take atoms that are present in 1kx5, and omit that are not present
                    self.alignment_sel_dict={}
                    al_sel=self.u.select_atoms(self.alignment_sel)
                    # !!!!!!!!!! ACHTUNG !!!!!!!!!!! #
                    # Здесь я поменял строки на то, чтобы сделать упорядоченную выборку!!!
                    # Так как выборка делается по сути по остаткам, то порядок должен быть правильным!
                    # строковая выборка по какому то не очень понятному алгоритму может получиться кривой
                    
                    #self.alignment_sel_dict['mobile']=' or '.join(['(index %d)'%a.index for a in al_sel.atoms if self.map_obj2ref_resids[a.segid].get((a.segid,a.resid),False)])
                    #self.alignment_sel_dict['reference']=' or '.join([f'(segid {self.map_obj2ref_resids[a.segid][(a.segid,a.resid)][0]} and resid {self.map_obj2ref_resids[a.segid][(a.segid,a.resid)][1]} and name {a.name})' for a in al_sel.atoms if self.map_obj2ref_resids[a.segid].get((a.segid,a.resid),False)])
                    self.alignment_sel_dict['mobile']=['(index %d)'%a.index for a in al_sel.atoms if self.map_obj2ref_resids[a.segid].get((a.segid,a.resid),False)]
                    self.alignment_sel_dict['reference']=[f'(segid {self.map_obj2ref_resids[a.segid][(a.segid,a.resid)][0]} and resid {self.map_obj2ref_resids[a.segid][(a.segid,a.resid)][1]} and name {a.name})' for a in al_sel.atoms if self.map_obj2ref_resids[a.segid].get((a.segid,a.resid),False)]
                if(ref=='first_frame'):
                    
                    refnuc=self.u
                    refnuc.trajectory[0]
                    self.alignment_sel_dict=self.alignment_sel
                alignment = align.AlignTraj(self.u, refnuc, select=self.alignment_sel_dict, in_memory=True,start=time[0],stop=time[1],step=time[2])
                alignment.run()      
                
                ############
                #Step 6. Identification of 3D mapping and components of DNA by structrual alignment - DNA dyad, etc. Also, set base pairing here
                ############
                #Dyad, иметь возможность строить графики параметров ДНК относительно диады.
                # иметь возможность select  определенные регионы ДНК относительно диады(?)
                # словарь dna_str_sel_dict={'DNAPOS0':'segid I J resid 0', 'DNAPOS1':'segid I and resid 1 or segid J and resid -1'
                # словарь маэпинга на 1kx5 map_ref2obj_dna3d={(segid,resid):(segid,resid)}
                   
                #Let's identify dyad
                #Sine nucleosome is aligned to NRF, dyad is at x~0, y>0
                # self.dyad_top=('I',0) #resid of dyad on top DNA strand
                # self.dyad_bot=('J',0) #resid of dyad on bottom DNA strand
                ##################################################################
                #!!!! TODO check for complementarity and adjust if nessesary !!!!#
                # test structures 6UPK, 1P3B, 5HQ2, 6MUO                         #
                # also  issue if there are more than 1 top and bot strand        #
                ##################################################################
                for k,v in self.components.items():
                    if v['entity']=='DNA':
                        if v['type']=='DNAtop':
                            #check if chain is close to dyad location
                            near_dyad_sel=self.u.select_atoms('segid %s and name C1\' and (prop abs x <= 15.0) and (prop y <= 55.0) and (prop y >= 25.0) and (prop abs z <= 15.0)'%k)
                            if len(near_dyad_sel)==0:
                                self.components[k]['type']=='DNAother'
                                logger.warning('Outlier DNAtop chain {k} found, dropping')
                            else:
                                s=self.u.select_atoms('segid %s and name C1\''%k)
                                #changed score formula to add a lot if atom is far from dyad x<15 25<y<55 z<15
                                #score= np.abs(s.positions[:,0])-np.sign(s.positions[:,1])*100
                                score= np.abs(s.positions[:,0]) + \
                                (200*((np.abs(s.positions[:,0])>15) | (25>s.positions[:,1]) | (s.positions[:,1]>55) | (np.abs(s.positions[:,2])>15)))



                                dist_from_center=np.sqrt(np.sum(s.positions**2,axis=1))
    #                             print(score)
                                self.dyad_top=(k,s.atoms[np.argmin(score)].resid)
                                # thats a cheaty way to get around altlocs missting stuf ETC
                                self.DNA_top_length=len(self.seqs[k]['strseq'])
                                self.DNA_left_length=self.dyad_top[1]-s.resids[0]
                                self.DNA_right_length=s.resids[-1]-self.dyad_top[1]
                                self.DNA_top_segid=k
                        elif v['type']=='DNAbot':
                            #check if chain is close to dyad location
                            near_dyad_sel=self.u.select_atoms('segid %s and name C1\' and (prop abs x <= 15.0) and (prop y <= 55.0) and (prop y >= 25.0) and (prop abs z <= 15.0)'%k)
                            if len(near_dyad_sel)==0:
                                self.components[k]['type']=='DNAother'
                                logger.warning('Outlier DNAbot chain {k} found, dropping')
                            else:
                                s=self.u.select_atoms('segid %s and name C1\''%k)
                                score= np.abs(s.positions[:,0])  + \
                                (200*((np.abs(s.positions[:,0])>15) | (25>s.positions[:,1]) | (s.positions[:,1]>55) | (np.abs(s.positions[:,2])>15)))
                                self.dyad_bot=(k,s.atoms[np.argmin(score)].resid)
                                self.DNA_bot_length=len(self.seqs[k]['strseq'])
                                self.DNA_bot_segid=k
                # that code checks if +-5 of dyad is complementare
                if check_dyad:
                    def fix_bot_dyad(dyad_top,dyad_bot,u):
                        def getsafe(array,index):
                            try:
                                return(array[index])
                            except IndexError:
                                return(None)

                        def check_compl(list1,list2):
                            compl_dict={'DA':'DT','DT':'DA','DG':'DC','DC':'DG',
                                       'ADE':'TYM','TYM':'ADE','GUA':'CYT','CYT':'GUA',
                                       'A':'T','T':'A','G':'C','C':'G',None:None}
                            return sum([True if ((base1 is None) or (base2 is None)) else base1==compl_dict.get(base2,None) for (base1,base2) in zip(list1,list2)])


                        topchain,topid=dyad_top
                        botchain,botid=dyad_bot
                        top_agrp=[u.select_atoms('segid %s and resnum %d'%(topchain,rid)) for rid in range(topid-5,topid+6)]
                        bot_agrp=[u.select_atoms('segid %s and resnum %d'%(botchain,rid)) for rid in range(botid+5,botid-6,-1)]
                        scores=[]
                        shifts=list(range(-4,5))
                        for shift in shifts:
                            if shift<0:
                                topsel=slice(None,shift)
                                botsel=slice(-shift,None)
                            elif shift>0:
                                topsel=slice(shift,None)
                                botsel=slice(None,-shift)
                            else:
                                topsel=slice(None,None)
                                botsel=slice(None,None)
                                #residues_3to1[
                            top_list=[getsafe(agrp.residues.resnames,0) for agrp in top_agrp[topsel]]
                            bot_list=[getsafe(agrp.residues.resnames,0) for agrp in bot_agrp[botsel]]

                            scores.append(check_compl(top_list,bot_list))


                        shift=shifts[scores.index(max(scores))]
                        return((dyad_bot[0],dyad_bot[1]+shift))
                    new_bot=fix_bot_dyad(self.dyad_top,self.dyad_bot,self.u)
                    self.dyad_bot= self.dyad_bot if new_bot is None else new_bot
                    
                #MAJOR CHANGES - idea is that we use bp_dict object, which is derived at structure initialisation
                # using allignment with rev compl bot strand
                self.bp_dict=get_base_pairs_dict(self)

    
    def write(self,path,step=1):
        """
        Write new pdb and coordinates 
        Currently what happens with the first frame is ambigous
        Path should be without extension, pdb and xtc will be appended.
        """
        sel=self.u.select_atoms('all')
        sel.write(path+'.pdb')
        self.u.trajectory[0]
        with mda.Writer(path+".xtc", sel.n_atoms) as W:
            for ts in tqdm(self.u.trajectory[self.time[0]:self.time[1]:self.time[2]*step]):
                W.write(sel)
        self.u.trajectory[0]
    
    def vmd_lv(self,write=True,show_1KX5=False,show_ref=False,step=1):
        """
        VMD local view code - will generate a code to download traj to your local computer and view with VMD
        Will also dump traj first.
        Experimental
        """
        if self.name is None:
            self.tmp_name=uuid.uuid4()
        else:
            self.tmp_name=self.name
        if not os.path.exists('tmp'):
            os.mkdir('tmp')
        if write:
            self.write('tmp/%s'%self.tmp_name,step=step)
        
        if show_1KX5 or show_ref=='1KX5':
            ref='1KX5_NRF.pdb'
        elif show_ref=='3LZ0':
            ref='3LZ0_NRF.pdb'
        else:
            ref=f"{self.tmp_name}.pdb"
        print("mkdir -p %s"%self.tmp_name)
        print("cd %s"%self.tmp_name)
        hn=socket.gethostname()
        DATA_PATH = pkg_resources.resource_filename('pynucl', 'data/')
        self.vmd_scr_type=''
#         print(DATA_PATH)
        print(f"rsync -avz {hn}:{DATA_PATH}/VMD_scripts/* .")
        print(f"rsync -avz {hn}:{os.path.abspath(f'tmp/{self.tmp_name}*')} .")
        print(f"vmd -e view_nucl{self.vmd_scr_type}.tcl -args {self.tmp_name}.pdb {self.tmp_name}.xtc {self.tmp_name} 0 0 1 1 1 1 1 0 0 {ref} \n\n")
        
        
        
    def view(self,**kwargs):
        '''
        Function shows preview of the nucleosome (chain names alike 1kx5) via NGLview and MDAnalysis
        args   - MDA universe or atom selection, or anything that can be parsed via mda.Universe()
        legend - show legend
        gui    - shows standard nglview gui (very limited)
        onlyNcp - show only ncp
        chconv - Dictionary with chain name conformity to map to 1kx5
        selection - mda selection string
        color  - color preset (bright or dull, or dictionary like color={'H3':"#94b4d1",'H4':"#94d19c",'H2A':"#d6d989",'H2B':"#d98989",'DNA':"#d6d6d6"})
        '''
        return view_nucl(self.u,**kwargs)
    
    def nucl_sel_expand(self,sel):
        return nucl_sel_expand(sel,self.nucl_elements)
            
class nucltrj(nuclstr):
    """
    Extends nuclstr to import trajectories
    """
    def __init__(self, topol, trj=None, name=None,topol_as_first_frame=True, **kwargs):
        if name is None:
            self.name=os.path.splitext(os.path.basename(topol))[0]
        if(isinstance(topol,mda.Universe)):
            super().__init__(topol, name, **kwargs)
        else:
            if(topol_as_first_frame and topol[-3:]=='pdb'):
                opened_trj=mda.Universe(topol,topol,trj)#,in_memory=True
            else:
                opened_trj=mda.Universe(topol,trj)#,in_memory=True
            super().__init__(opened_trj, name, **kwargs)
            
            
    

    

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

