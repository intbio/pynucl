# -*- coding: utf-8 -*-

"""
#This module provides routines to identify the entities of nolecule elements (histone, dna e.t.c.).



Copyright 2020, Alexey Shaytan, Grigoriy Armeev, Moscow Sate Univeristy

"""

#from .base import *
#from nucl_ref_frame_aln import nucl_align
#from .nucl_meta_select import create_elem_dict,nucl_sel_expand
#from .view_nucl import view_nucl
#from .seq_utils import *
#from .hist_features import *

import MDAnalysis as mda
from MDAnalysis.analysis import align
import pandas as pd
import numpy as np
#import tempfile
#from io import StringIO
#import os
#from tqdm.auto import tqdm
#import uuid
#import socket
#import pkg_resources

# plan for 06.05.2020
# - create a mockup for all histone-related stuff -done 
# - create a database file from histonbv e.t.c with all standard histone sequences -done 08.05
#plan for 08.05
# - create a database for dna variants
# - write logick for dna detection

# the goal
'''
{'I':{'entity':'DNA','type':'DNAtop','variant':'alphasat','side':0},
 'J':{'entity':'DNA','type':'DNAbot','variant':'alphasat','side':0},
 'A':{'entity':'histone','type':'H3','variant':'canonicalH3','side':1},
 'E':{'entity':'histone','type':'H3','variant':'canonicalH3','side':2},
 'B':{'entity':'histone','type':'H4','variant':'canonicalH4','side':1},
 'F':{'entity':'histone','type':'H4','variant':'canonicalH4','side':2},
 'C':{'entity':'histone','type':'H2A','variant':'canonicalH2A','side':1},
 'G':{'entity':'histone','type':'H2A','variant':'canonicalH2A','side':2},
 'D':{'entity':'histone','type':'H2B','variant':'canonicalH2B','side':1},
 'H':{'entity':'histone','type':'H2B','variant':'canonicalH2B','side':2}}'''



import logging
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)s %(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

#what do we have:
#struct,name=None, format=None, time=(0,None,1), ref='1KX5_NRF', fullseqs=None, debug=False,skipaln=False
#struct - an MDanalysis universe (single structrue or trajectory)
#fullseqs - information about the full length sequences of the corresponding protein and DNA chains in the structure.
#        can be: 'PDB_ID' (e.g. 1KX5) - will try to get info directly from PDB fasta header files (SEQREC records), this might be or might not be the full sequence


from .utils.structure_constants import WaterNames,IonNames,DnaBases,RnaBases,NucleicBackboneAtoms
from .hist_ss import identify_hist_type
from .hist_features import hist_ident
from .seq_utils import seqrec_from_pdb,seq_from_mda,get_overhangL,get_overhangR
from .nucl_meta_select import create_elem_dict,nucl_sel_expand


def detect_entities(struct=None,fullseqs=None,blast_evalue=100,blast_positive_ratio=0.8):
    #using set to remove duplicates
    entity_dict={segid:{} for segid in sorted(set(struct.segments.segids))}
    
    for segid,segid_dict in entity_dict.items():
        segid_sel  =struct.select_atoms(f'segid {segid}')
        nucleic_sel=segid_sel.select_atoms('nucleic')
        protein_sel=segid_sel.select_atoms('protein')
        water_sel  =segid_sel.select_atoms(f'resname {" ".join(WaterNames)}')
        ion_sel    =segid_sel.select_atoms(f'resname {" ".join(IonNames)}')
        other_sel  =segid_sel - (nucleic_sel + protein_sel + water_sel + ion_sel)


        if (len(protein_sel) > 0):
            # here we have to decide if that thing is histone or other thing
            seq=seq_from_mda(protein_sel,segid)
            hist_type,hist_variant,hist_oganism=identify_hist_type_variant_organism(seq,
                                                                                    evalue=blast_evalue,
                                                                                    positive_ratio=blast_positive_ratio)
            if hist_type is None:
                segid_dict['entity']='other_protein'
            else:
                segid_dict['entity']='histone'
                segid_dict['type']=hist_type
                segid_dict['variant']=hist_variant
                segid_dict['organism']=hist_oganism
            segid_dict['sequence']=seq

        elif (len(nucleic_sel) > 0) and (len(protein_sel) == 0):
            # here we have to detect if that is DNA or RNA, that may be tricky
            # lets try to detect if it is RNA by base name or oxygen atom
            has_rna_names   =bool(len(nucleic_sel.select_atoms(f'resname {" ".join(RnaBases)}')))

            rna_specific_atoms=nucleic_sel.select_atoms('name O2*')
            # mdanalysis understands * as wildcard, so have to check
            has_rna_specific_atoms=any([True for name in rna_specific_atoms.names if name in ["O2'","O2*"]])
            
            if has_rna_names and has_rna_specific_atoms:
                segid_dict['entity']='RNA'
            else:
                segid_dict['entity']='DNA'
                seq=seq_from_mda(nucleic_sel,segid)
                variant,type=identify_DNA_chain(seq,evalue=blast_evalue)
                segid_dict['type']=type
                segid_dict['variant']=variant
                segid_dict['sequence']=seq


            # top and bot strabd deduction idea - by contacts with h3 h4
        #elif (len(nucleic_sel) > 0) and (len(protein_sel) > 0):
        #    segid_dict['entity']='mixed'
        else:
            segid_dict['entity']='unknown'
        segid_dict['has_water']=bool(len(water_sel))
        segid_dict['has_ions']=bool(len(ion_sel))
        segid_dict['has_other']=bool(len(other_sel))
        if segid_dict['entity']=='histone':
            # that one is needed as we will determine the side later
            segid_dict['side']=1
        else:
            segid_dict['side']=0
    entity_dict=check_dna_top_and_bottom_strands(entity_dict)
    return(entity_dict)

import itertools
from Bio import pairwise2
import itertools
from Bio import pairwise2
def check_dna_top_and_bottom_strands(entity_dict):
    def get_chains_complementarity(seq1,seq2):
        alignments = pairwise2.align.globalxx(Seq(str(seq1)),Seq(str(seq2)).reverse_complement(),penalize_end_gaps=False)
        return(2*alignments[0][-1]/(len(seq1)+len(seq2)))
    
    dna_chain_list=[key for key,item in entity_dict.items() if item['entity']=='DNA']
    if len(dna_chain_list)!=0:
        top_segid=bot_segid=None
        #check that there is only singular top and bottom strand
        top_strands=[entity_dict[segid]['type']=='DNAtop' for segid in dna_chain_list]
        bot_strands=[entity_dict[segid]['type']=='DNAbot' for segid in dna_chain_list]
        if (sum(top_strands)+sum(bot_strands))==2:
            if (sum(top_strands)==1) and (sum(bot_strands)==1):
                top_strands.index(True)
                top_segid=dna_chain_list[top_strands.index(True)]
                bot_segid=dna_chain_list[bot_strands.index(True)]
                logger.debug('Top and bottom strands detected')

            elif (sum(top_strands)==2) and (sum(bot_strands)==0):
                top_segid,bot_segid=[dna_chain_list[i] for i in range(len(dna_chain_list)) if top_strands[i]]            
                logger.debug(f'Two top strands detected, likely palindromic sequence, reassigning by chain name Chains {top_segid} and {bot_segid}')

            elif (sum(top_strands)==0) and (sum(bot_strands)==2):
                top_segid,bot_segid=[dna_chain_list[i] for i in range(len(dna_chain_list)) if bot_strands[i]]
                logger.debug(f'Two bot strands detected, likely palindromic sequence, reassigning by chain name Chains {top_segid} and {bot_segid}')

        elif (sum(top_strands)==1) and (sum(bot_strands)==0):
            logger.debug('Single top strand detected, finding complementary bot')
            top_segid=dna_chain_list[top_strands.index(True)]
            dna_chain_list.pop(dna_chain_list.index(top_segid))
            score=[]
            for segid2 in dna_chain_list:
                score.append(get_chains_complementarity(entity_dict[top_segid]['sequence'],entity_dict[segid2]['sequence']))
            if max(score)>0.85:
                bot_segid=dna_chain_list[score.index(max(score))]
                logger.debug(f'Chains {top_segid} and {bot_segid} are complementary assigning')


        elif (sum(top_strands)==0) and (sum(bot_strands)==1):
            logger.debug('Single bot strand detected, finding complementary top')
            bot_segid=dna_chain_list[bot_strands.index(True)]
            dna_chain_list.pop(dna_chain_list.index(bot_segid))
            score=[]
            for segid1 in dna_chain_list:
                score.append(get_chains_complementarity(entity_dict[segid1]['sequence'],entity_dict[bot_segid]['sequence']))
            if max(score)>0.85:
                top_segid=dna_chain_list[score.index(max(score))]
                logger.debug(f'Chains {top_segid} and {bot_segid} are complementary assigning')


        else:
            logger.warning('No top and bot strands detected from database, finding two complementary chains and assigning by chain name')
            dna_chain_combinations=list(itertools.combinations(dna_chain_list,2))
            score=[]
            for segid1,segid2 in dna_chain_combinations:
                score.append(get_chains_complementarity(entity_dict[segid1]['sequence'],entity_dict[segid2]['sequence']))
            if max(score)>0.85:
                top_segid,bot_segid=dna_chain_combinations[score.index(max(score))]
                logger.warning(f'Chains {top_segid} and {bot_segid} are complementary assigning top and bot by chain name')
            for segid1,segid2 in dna_chain_combinations:
                if (segid1=='I') and (segid2=='J'):
                    top_segid='I'
                    bot_segid='J'
                    logger.warning(f'Found I and J, assigning blindly')
                    break
        if not(top_segid is None) and not (bot_segid is None):
            entity_dict[top_segid]['type']='DNAtop'
            entity_dict[bot_segid]['type']='DNAbot'
        else:
            logger.warning('Could not find chains with significant complementarity')
    # TODO check that the sewuence detected is the same )
    return (entity_dict)
    
    

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform,cdist
import seaborn as sns
from itertools import combinations,product,chain
from pymolint import mol_int
import copy
def check_nucleosomes_sides(entity_dict,elements_dict,struct,seq_features,pdbid=None):
    components=entity_dict
    u=struct
    coms={}
    hist_folds=u.select_atoms(nucl_sel_expand('hist_folds',elements_dict))
    for k in components.keys():
        if ( components[k]['entity']=='histone') and (components[k]['type']!='H1') and (len(seq_features[k])>1):
            coms[k]=hist_folds.select_atoms(f'name CA and segid {k}').center_of_geometry()

    com_list=[item for key,item in coms.items()]
    segid_list=[key for key,item in coms.items()]

    pairthreshold=10 #A
    H3H4_pairs=[]
    H2AH2B_pairs=[]
    other_pairs=[]
    for pair,distance in zip(combinations(segid_list,2),pdist(com_list)):
        if distance<=pairthreshold:
            hist1,hist2=components[pair[0]]['type'],components[pair[1]]['type']
            if (set((hist1,hist2)).issubset(['H3','H4']) and (hist1!=hist2)):      
                H3H4_pairs.append(pair)
            elif (set((hist1,hist2)).issubset(['H2A','H2B']) and (hist1!=hist2)):
                H2AH2B_pairs.append(pair)
            else:
                other_pairs.append(pair)
    coms={}
    alpha3=u.select_atoms(nucl_sel_expand('H3 and alpha3',elements_dict))
    for pair in H3H4_pairs:
        coms[pair]=alpha3.select_atoms(f'name CA and segid {" ".join(pair)}').center_of_geometry()
    com_list=[item for key,item in coms.items()]
    pair_list=[key for key,item in coms.items()]

    H3H3threshold=16 #A
    H3H4_H3H4_list=[]
    for pair,distance in zip(combinations(pair_list,2),pdist(com_list)):
        if distance<=H3H3threshold:
            H3H4_H3H4_list.append(list(pair))
    H3H4_H3H4_list    
    #find H2AH2B_H2AH2B tetrameres
    coms={}
    loopL1=u.select_atoms(nucl_sel_expand('H2A and loopL1',elements_dict))
    for pair in H2AH2B_pairs:
        coms[pair]=loopL1.select_atoms(f'name CA and segid {" ".join(pair)}').center_of_geometry()

        com_list=[item for key,item in coms.items()]
    pair_list=[key for key,item in coms.items()]

    H2AH2Bthreshold=20 #A
    H2AH2B_H2AH2B_list=[]
    for pair,distance in zip(combinations(pair_list,2),pdist(com_list)):
        if distance<=H2AH2Bthreshold:
            H2AH2B_H2AH2B_list.append(list(pair))
    H2AH2B_H2AH2B_list
    
    # find hemisomes
    H3H4_H2AH2B_same_side_list=[]
    H3H4_H2AH2B_opposed_side_list=[]
    H3H4_H2AH2B_outlier_list=[]
    beta3_H4=u.select_atoms(nucl_sel_expand('name CA and beta3 and H4',elements_dict))
    beta3_H2A=u.select_atoms(nucl_sel_expand('name CA and beta3 and H2A',elements_dict))

    alpha3_H4=u.select_atoms(nucl_sel_expand('name CA and alpha3 and H4',elements_dict))
    alpha3_H2B=u.select_atoms(nucl_sel_expand('name CA and alpha3 and H2B',elements_dict))

    opposed_side_threshold=10
    same_side_threshold=20
    for h3h4,h2ah2b in product(H3H4_pairs,H2AH2B_pairs):
        #opposed (H4 or H2A) and beta3
        beta3_H4_com = beta3_H4.select_atoms(f'segid {" ".join(  h3h4)}').center_of_geometry()
        beta3_H2A_com=beta3_H2A.select_atoms(f'segid {" ".join(h2ah2b)}').center_of_geometry()
        opposed_side_dist=np.linalg.norm(beta3_H4_com-beta3_H2A_com)

        #same (H4 or H2B) and alpha3
        alpha3_H4_com = alpha3_H4.select_atoms(f'segid {" ".join(  h3h4)}').center_of_geometry()
        alpha3_H2B_com=alpha3_H2B.select_atoms(f'segid {" ".join(h2ah2b)}').center_of_geometry()

        same_sede_dist=np.linalg.norm(alpha3_H4_com-alpha3_H2B_com)
        #same_sede_dist==0 is needed as selection can be empty
        if (opposed_side_dist<=opposed_side_threshold) and ((same_sede_dist >= same_side_threshold) or(same_sede_dist==0) ):
            H3H4_H2AH2B_opposed_side_list.append([h3h4,h2ah2b])
        elif ((opposed_side_dist>=opposed_side_threshold) or (opposed_side_dist==0)) and (same_sede_dist <= same_side_threshold):
            H3H4_H2AH2B_same_side_list.append([h3h4,h2ah2b])
        elif (opposed_side_dist<=opposed_side_threshold) and (same_sede_dist <= same_side_threshold):
            H3H4_H2AH2B_outlier_list.append([h3h4,h2ah2b])
    if len(H3H4_H2AH2B_outlier_list)!=0:
        print('Warning, detected strange geometry during hemisome detection ',pdbid)

    #we will have to order it after
    nlp_list=[]
    #1. add all tetrameres
    nlp_list.extend(copy.deepcopy(H3H4_H3H4_list))
    #2. add loose h3h4 dimers
    for H3H4_pair in H3H4_pairs:
        if not(H3H4_pair in chain.from_iterable(nlp_list)):
            nlp_list.append([H3H4_pair]) 
    #3. add same side dimers

    for H3H4_H2AH2B in H3H4_H2AH2B_same_side_list:
        for nlp in nlp_list:
            if (H3H4_H2AH2B[0] in nlp):   
                nlp.append(H3H4_H2AH2B[1])
            elif (H3H4_H2AH2B[1] in nlp):
                nlp.append(H3H4_H2AH2B[0])

    # #4. check that opposed dimers are in nlps
    for H3H4_H2AH2B in H3H4_H2AH2B_opposed_side_list:
        if not set(H3H4_H2AH2B).issubset(chain.from_iterable(nlp_list)):
            print('Loose H3H4 H2AH2B tetramer detected segments:',H3H4_H2AH2B,pdbid)
            for nlp in nlp_list:
                if (H3H4_H2AH2B[0] in nlp):   
                    nlp.append(H3H4_H2AH2B[1])
                elif (H3H4_H2AH2B[1] in nlp):
                    nlp.append(H3H4_H2AH2B[0])
    # #5. check that h2a h2b dimers are in nlps
    for H2AH2B_H2AH2B in H2AH2B_H2AH2B_list:
        if not set(H2AH2B_H2AH2B).issubset(chain.from_iterable(nlp_list)):
            print('Loose H3H4 H2AH2B tetramer detected segments:',H2AH2B_H2AH2B,pdbid)

    # #6. create dictionary of nlps with anotations
    # nlp_list

            
    nlp_dict={}
    for i,nlp in enumerate(nlp_list):
        indict={}
        indict['all']=list(chain.from_iterable(nlp))
        indict['all'].sort()
        indict['H3']=[segid for segid in indict['all'] if components[segid]['type']=='H3']
        indict['H4']=[segid for segid in indict['all'] if components[segid]['type']=='H4']
        indict['H2A']=[segid for segid in indict['all'] if components[segid]['type']=='H2A']
        indict['H2B']=[segid for segid in indict['all'] if components[segid]['type']=='H2B']
        indict['H3_H4']=[dimer for dimer in nlp if dimer in H3H4_pairs]
        indict['H2A_H2B']=[dimer for dimer in nlp if dimer in H2AH2B_pairs]
        indict['sides']=[[dimer] for dimer in indict['H3_H4']]
        for dimer in indict['sides']:
            for same_side_tetramer in H3H4_H2AH2B_same_side_list:
                if dimer[0] in same_side_tetramer:
                    dimer.extend([other_dimer for other_dimer in same_side_tetramer if other_dimer not in dimer])

        if len(indict['sides'])==1 and len(H3H4_H2AH2B_opposed_side_list)==1:
            if indict['sides'][0][0] in H3H4_H2AH2B_opposed_side_list[0]:
                indict['sides'].append([other_dimer for other_dimer in H3H4_H2AH2B_opposed_side_list[0] if other_dimer not in dimer])

        if len(indict['sides'])==2:
            if len(indict['H2A_H2B'])==2:
                indict['type']='octamer'
            elif len(indict['H2A_H2B'])==1:
                if len(indict['H3_H4'])==2:
                    indict['type']='hexamer'
                elif len(indict['H3_H4'])==1:
                    indict['type']='split_octasome'
            elif len(indict['H2A_H2B'])==0:
                indict['type']='H3_H4_tetramer'
        elif len(indict['sides'])==1:
            if len(indict['sides'][0])==2:
                indict['type']='hemisome'
            elif len(indict['H3_H4'])==1:
                indict['type']='H3_H4_dimer'
            elif len(indict['H2A_H2B'])==1:
                indict['type']='H2A_H2B_dimer'
        nlp_dict[i]=indict
    nlp_dict   

    # (component['entity']=='DNA') and (component['type']=='DNAtop' ) order is essential, as there are components without type kwrd
    top_dna_segid= [segid for segid,component in components.items() if (component['entity']=='DNA') and (component['type']=='DNAtop' )][0]
    #print(top_dna_segid)
    for nlpid,nlp in nlp_dict.items():
        nlp['sideids']=[]
        for sideid,side_sigids in enumerate(nlp['sides'],1):
            side_segid_list=list(chain.from_iterable(side_sigids))
            h3_one_sel=u.select_atoms(f'segid {" ".join(side_segid_list)}')

            point1=h3_one_sel.select_atoms(nucl_sel_expand('name CA and H3 and loopL1',elements_dict)).center_of_geometry()
            point2=h3_one_sel.select_atoms(nucl_sel_expand('name CA and H3 and mgarg1',elements_dict)).center_of_geometry()
            point3=h3_one_sel.select_atoms(nucl_sel_expand('name CA and H3 and loopL2',elements_dict)).center_of_geometry()


            top_chanin_sel=u.select_atoms(f'segid {top_dna_segid} and name P')
            top_chanin_sel_pos=top_chanin_sel.positions
            top_chanin_sel_resids=top_chanin_sel.resids
            points=[point for point in [point1,point2,point3] if len(point)==3]
            if len(points)==2:
                print('Warning either some histone features missing or mapping is wrong')
            if len(points)>=2:
                dist_mat=cdist(points,top_chanin_sel_pos)
                key_rid=top_chanin_sel_resids[np.argmin(dist_mat,axis=1)]
                if np.all(np.diff(key_rid)>=0):
                    nlp['sideids'].append(1)
                else:
                    nlp['sideids'].append(2)
            else:
                if len(side_segid_list)==2 and set([components[segid]['type'] for segid in side_segid_list]).issubset(['H2A','H2B']):
                    nlp['sideids'].append(1 if nlp['sideids'][-1]==2 else 2)
                    print('guessed side for chains '+" ".join(side_segid_list))
                else:
                    print('cant detect side for chains '+" ".join(side_segid_list))

    # #6. create dictionary of nlps with anotations
    # nlp_list
    nlp_list
    for nlpid, nlp in nlp_dict.items():
        for side_sigids,sideid in zip(nlp['sides'],nlp['sideids']):
            for segid in chain.from_iterable(side_sigids):
                component=components[segid]
                component['side']=sideid
                component['nlp_id']=nlpid
                component['nlp_type']=nlp['type']
    return components

from .hist_features import hist_features
def detect_entities_independent(struct=None,fullseqs=None,blast_evalue=100,blast_positive_ratio=0.8):
    components=detect_entities(struct=struct,fullseqs=fullseqs,blast_evalue=blast_evalue,blast_positive_ratio=blast_positive_ratio)
    seqs={}
    seqs={k:{} for k in components.keys()}#populate with keys
    u=struct
    selection=u.select_atoms(f'(protein or nucleic) and (not resname {" ".join(IonNames)})')
    for k in components.keys():
        if(components[k]['entity']=='histone'):
            segment=selection.select_atoms(f'segid {k}')
            strseq=seq_from_mda(segment,k)
            resids=segment.residues.resids
            seqs[k]['strseq']=strseq[0]+''.join(['X'*(diff-1) + strseq[i+1] for i,diff in enumerate(np.diff(resids))])
    ##Step 2.2. get fullsequences
    if (fullseqs is None): #inherit strseqs
        for k in components.keys():
            if(components[k]['entity']=='histone'):
                #self.seqs[k]['fullseq']=self.seqs[k]['strseq']
                seqs[k]['fullseq']=seqs[k]['strseq']
                seqs[k]['overhangL']=seqs[k]['overhangR']=0
                seqs[k]['resid_start']=u.select_atoms('segid %s and (protein or nucleic)'%k).residues.resids[0]
                components[k]['sequence']=seqs[k]['fullseq']
                components[k]['resid_start']=int(seqs[k]['resid_start'])
    else:
        if isinstance(fullseqs, str) and len(fullseqs)==4:
            for k in components.keys():
                if(components[k]['entity']=='histone'):
                    seqs[k]['fullseq']=seqrec_from_pdb(fullseqs,k)
        elif isinstance(fullseqs, dict):
            for k in components.keys():
                if(components[k]['entity']=='histone'):
                    seqs[k]['fullseq']=fullseqs[k]


        for k in self.components.keys():
            if(components[k]['entity']=='histone'):
                seqs[k]['overhangL']=get_overhangL(seqs[k]['fullseq'],seqs[k]['strseq'])
                seqs[k]['overhangR']=get_overhangR(seqs[k]['fullseq'],seqs[k]['strseq'])
                seqs[k]['resid_start']=u.select_atoms('segid %s and (protein or nucleic)'%k).residues.resids[0]
                components[k]['sequence']=seqs[k]['fullseq']
                components[k]['resid_start']=int(seqs[k]['resid_start'])

    seq_features={}
    for k in components.keys():
        if(components[k]['entity']=='histone'):
            if components[k]['type']!='H1':
                seq_features[k]=hist_features(seqs[k]['fullseq'])
    elements_dict=create_elem_dict(components,seqs,seq_features)
    components=check_nucleosomes_sides(components,elements_dict,u,seq_features)
    return(components)

### Helper - should be transfered somewhere sequence-related
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline,NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import uuid,os,re

from .utils.memory_tempfile import MemoryTempfile
tempfile = MemoryTempfile()
try:
    this_dir=os.path.split(__file__)[0]
    logger.debug(this_dir)
    DATA_PATH = os.path.join(this_dir, "data")
except:
    DATA_PATH=None

def identify_hist_type_variant_organism(test_seq,evalue=100,positive_ratio=0.8):
    '''
    performs local blast search via histonedb seeds database
    input:
    test_seq - either BIO.Seq object, or string
    evalue - blast evalue threshold (default 100)
    positive_ratio - N positive matches to queue sequence length (default 0.8)
    returns:
    str(hist_type)
    str(hist_variant)
    str(hist_oganism)
    #Done - use tempfile like in 3dna wrapper
    #TODO - may be better to use pairwise against prealligned sequences
    '''
    # Merge all sequences to one
    histone_seeds={'H3_fasta':os.path.join(DATA_PATH,'histonedb_seq_seeds','H3.fasta'),
                   'H4_fasta':os.path.join(DATA_PATH,'histonedb_seq_seeds','H4.fasta'),
                   'H2A_fasta':os.path.join(DATA_PATH,'histonedb_seq_seeds','H2A.fasta'),
                   'H2B_fasta':os.path.join(DATA_PATH,'histonedb_seq_seeds','H2B.fasta'),
                   'H1_fasta':os.path.join(DATA_PATH,'histonedb_seq_seeds','H1.fasta')}
    histone_seq_records=[]
    for histone_seed in histone_seeds.values():
        histone_seq_records.extend(list(SeqIO.parse(histone_seed, "fasta")))
    
    #remove gaps from MSA
    for record in histone_seq_records:
        record.seq=record.seq.ungap("-")
    
    n1=str(uuid.uuid4())
    n2=str(uuid.uuid4())
    seqLen=len(test_seq)
    if (not isinstance(test_seq,Seq)) and (isinstance(test_seq,str)):
        test_seq=Seq(test_seq)
    elif (not isinstance(test_seq,Seq)): raise TypeError("Test sequence must be either Bio.Seq or string")
    
    with tempfile.TemporaryDirectory() as TEMP:
        SeqIO.write([SeqRecord(test_seq,id='Query',name='Query')],os.path.join(TEMP,n2+'.fasta'),'fasta')

        SeqIO.write(histone_seq_records, os.path.join(TEMP,n1+'.faa'), "fasta")
        os.system('makeblastdb -dbtype prot -in %s.faa -out %s.db > /dev/null'%(os.path.join(TEMP,n1),
                                                                                os.path.join(TEMP,n1)))


        blastp_cline = NcbiblastpCommandline(query=os.path.join(TEMP,n2+'.fasta'),
                                             db=os.path.join(TEMP,n1+'.db'),
                                             evalue=evalue,outfmt=5,
                                             out=os.path.join(TEMP,n1+'.xml'))
        stdout, stderr = blastp_cline(cwd=TEMP)

        blast_record = NCBIXML.read(open(os.path.join(TEMP,n1+'.xml'),'r'))

        sname=list()
        evalue=list()
        score=list()
        hsp_list=list()
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:           
                sname.append(alignment.title)
                evalue.append(hsp.expect)
                # this is not blast score that is persenatage of identity of the queue seq to database
                score.append(hsp.positives/seqLen ) 
                hsp_list.append(hsp)
                # length_list.append(alignment.length)
        hist_variant=hist_type=hist_oganism=None
        if len(evalue)>5: # there should be a lot of blast hists, lrts define threshold as 5 (no reason)
            min_index=evalue.index(min(evalue))

            if score[min_index]>positive_ratio: # 80% of the queue sequence should be alligned
                hist_identified=sname[min_index].split()[1]
                hsp=hsp_list[min_index]
                hist_detected=hist_identified.split('|')
                hist_variant=hist_detected[0]
                hist_type_re=re.search('H1|H2A|H2B|H3|H4',hist_variant)
                hist_type=hist_variant[slice(*hist_type_re.span())]
                hist_oganism=hist_detected[1]



        # length=length_list[evalue.index(min(evalue))]
        #os.system("rm %s.faa %s.db.phr %s.db.pin %s.db.psq %s.fasta %s.xml"%(n1,n1,n1,n1,n2,n1))
    return hist_type,hist_variant,hist_oganism


def identify_DNA_chain(test_seq,evalue=100):
    '''
    performs local blast search via known NPS
    input:
    test_seq - either BIO.Seq object, or string
    evalue - blast evalue threshold (default 100)
    returns:
    str(sequence_name)
    str(direction) ('top', 'bot')
    # TODO: check for telomeric sequence is stupid
    '''
    # Merge all sequences to one
    nps_seeds={'all':os.path.join(DATA_PATH,'positioning_sequences','sequences.fasta')}
    db_seq_records=[]
    for nps_seed in nps_seeds.values():
        db_seq_records.extend(list(SeqIO.parse(nps_seed, "fasta")))
    nps_seq_records=[]
    for seq_rec in db_seq_records:
        # TERRIBLE SOLUTION REFACTOR!!!!
        #dyad_locations=seq_rec.description.split('|')[1].split()[1].split(',')
        nps_seq_records.append(seq_rec)
        #for dyad in dyad_locations:
        #    nps_seq_records.append(seq_rec)
            #nps_seq_records[-1].seq=seq_rec.seq[int(dyad)-70:int(dyad)+70]
            #nps_seq_records.append(nps_seq_records[-1].reverse_complement())
            #nps_seq_records[-1].id=seq_rec.id+'_rev_comp'
            #nps_seq_records[-1].name=seq_rec.name+'_rev_comp'
            #nps_seq_records[-1].description=seq_rec.description
            

    #remove gaps from MSA
    for record in nps_seq_records:
        record.seq=record.seq.ungap("-")
    
    n1=str(uuid.uuid4())
    n2=str(uuid.uuid4())
    if (not isinstance(test_seq,Seq)) and (isinstance(test_seq,str)):
        test_seq=Seq(test_seq)
    else: raise TypeError("Test sequence must be either Bio.Seq or string")
    with tempfile.TemporaryDirectory() as TEMP:
        SeqIO.write([SeqRecord(test_seq,id='Query',name='Query')],os.path.join(TEMP,n2+'.fasta'),'fasta')

        SeqIO.write(nps_seq_records, os.path.join(TEMP,n1+'.faa'), "fasta")
        os.system('makeblastdb -dbtype nucl -in %s.faa -out %s.db > /dev/null'%(os.path.join(TEMP,n1),
                                                                                os.path.join(TEMP,n1)))

        blastn_cline = NcbiblastnCommandline(query=os.path.join(TEMP,n2+'.fasta'),
                                             db=os.path.join(TEMP,n1+'.db'),
                                             evalue=evalue,outfmt=5,
                                             strand='both',word_size=20, perc_identity=90,
                                             out=os.path.join(TEMP,n1+'.xml'))
        stdout, stderr = blastn_cline(cwd=TEMP)

        blast_record = NCBIXML.read(open(os.path.join(TEMP,n1+'.xml'),'r'))

        sname=list()
        evalue=list()
        hsp_list=list()
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sname.append(alignment.title)
                evalue.append(hsp.expect)
                hsp_list.append(hsp)

                # length_list.append(alignment.length)
        if len(evalue)>0: 
            nps_identified=sname[evalue.index(min(evalue))].split()[1]
            strand=hsp_list[evalue.index(min(evalue))].strand[1]
            if strand=='Plus':
                direction='DNAtop'
            elif strand=='Minus':
                direction='DNAbot'
            
        elif 'TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG' in str(test_seq):
            direction='DNAtop'
            nps_identified='telomeric_human'
        elif 'CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA' in str(test_seq):
            direction='DNAbot'
            nps_identified='telomeric_human'
        else:
            direction=nps_identified=None

        #os.system("rm %s.faa %s.db.nhr %s.db.nin %s.db.nsq %s.fasta %s.xml"%(n1,n1,n1,n1,n2,n1))
    return nps_identified,direction


from Bio import Align
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
def get_base_pairs_dict(p,verbose=False):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -0.5
    aligner.open_gap_score = -2.0
    aligner.end_open_gap_score = 0
    dnatop=Seq(str(p.seqs[p.dyad_top[0]]['strseq']), generic_dna)
    dnabot=Seq(str(p.seqs[p.dyad_bot[0]]['strseq']), generic_dna)
    dnabot.reverse_complement()
    #alignments = aligner.align(dnatop, dnabot.reverse_complement())
    alignments = aligner.align(str(dnatop), str(dnabot.reverse_complement()))
    best_score=0
    best_aln='no'
    i=0
    for a in sorted(alignments):
        if(a.score>best_score):
            best_score=a.score
            best_aln=a
            i+=1
        if i>100:
            break
    if verbose:
        print(best_aln)
    top_resid_offset=p.seqs[p.dyad_top[0]]['resid_start']
    bot_resid_offset=p.seqs[p.dyad_bot[0]]['resid_start'] + p.DNA_bot_length -1
    basepairs={'top':{},'bot':{}}
    for top_strip,bot_strip in zip(*best_aln.aligned):
        for i_top,i_bot in zip(list(range(*top_strip)),list(range(*bot_strip))):
            if not 'X' in [best_aln.target[i_top],best_aln.query[i_bot]]:
                top_base=top_resid_offset+i_top
                bot_base=bot_resid_offset-i_bot
                basepairs['top'][top_base]=bot_base
                basepairs['bot'][bot_base]=top_base
    return basepairs


