# -*- coding: utf-8 -*-
# The meta selection routines for nucleosome structural elements
# What we want to be able to do:
# To our analysis routines we'd like to be able to pass meta selection rules
# i.e. select histone folds of hitsones H3 and H4 on the proximal side.
# The easiest language for that is something like:
# select (side) and (H3 or H4) and (hist_folds)
# We need to extend mdanalysis syntax with
# side1 H3 H4 hist_folds alpha1 etc definitions
# this means every nuclstr object should have a predefined dictionary of selections
# and their should be a parser to expand the language into MDAnalysis selections.


#Here we should have helper functions to pupulate the nuclstr.nucl_elements dictionary.


#And helper functions to expand the meta selection syntaxis 

from .hist_features import *
import pandas as pd


def create_elem_dict(components,seqs,seq_features):
    """
    Given:
    1) a dictionary of the univese components (from pynucl)
    2) their sequences (from pynucl - the seqs dictionary)
    3) histone elements features.
    return a dictionary of nucleosome elements via MDAnalysis selection syntax.
    """
    
    
    #Step 1. Generate a specific dictionary for 
    # secondary structure (ss) elements in every segid
    
    ss_seq_dict={}
    for k,v in components.items():
        ss_seq_dict[k]={}
        for i in seq_features.get(k,[]):
            ss_seq_dict[k][i.id]=[i.location.start,i.location.end-1]
        #ss_seq_dict[k]=ss4seq(seqs[k]['strseq'])
    
    
    #Step 2. populate nucl_elements with entity related elements - sides and histones (DNA?)
    nucl_elements={'side1':'(not all)','side2':'(not all)','H3':'(not all)','H4':'(not all)','H2A':'(not all)','H2B':'(not all)'}
    for k,v in components.items():
        if v['entity']=='histone':
            nucl_elements['side%d'%v['side']]=nucl_elements['side%d'%v['side']]+'or (segid %s)'%k
            nucl_elements[v['type']]=nucl_elements[v['type']]+'or (segid %s)'%k
    
    
    #Step 3. Populate with ss_elements
    for i in basic_elem_terms:
        nucl_elements[i]='(not all)'
        
    for k,v in components.items():
        if v['entity']=='histone':
            for k2,v2 in ss_seq_dict[k].items():
                nucl_elements[k2]=nucl_elements[k2]+'or (segid %s and resid %d:%d)'%(k,v2[0]+1,v2[1]+1)

    #Step 4. Populate with groups (super selections)
    for k,v in elem_groups.items():
        nucl_elements[k]='(not all)'
        for i in v:
            nucl_elements[k]=nucl_elements[k]+'or (%s)'%nucl_elements[i]
                
                
    #as a precausion wrap everything in parenthesis
    for k,v in nucl_elements.items():
        nucl_elements[k]="(%s)"%v
    
    
    return nucl_elements


def nucl_sel_expand(selection,elem_dict):
    """
    Given an element dict={'def':'resid 1:2'} and selection='name CA and def'
    return 'name CA and resid 1:2
    """
    sel=selection
    for k, v in elem_dict.items():
        sel = sel.replace(k, v)
    
    return sel
