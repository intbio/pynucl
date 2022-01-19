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
            if 'resids' in seqs[k].keys():
                #todo write exceptions!
                ss_seq_dict[k][i.id]=[seqs[k]['resids'][i.location.start],seqs[k]['resids'][i.location.end-1]]
            else:
                ss_seq_dict[k][i.id]=[seqs[k]['resid_start']-seqs[k]['overhangL']+i.location.start,seqs[k]['resid_start']-seqs[k]['overhangL']+i.location.end-1]
        #ss_seq_dict[k]=ss4seq(seqs[k]['strseq'])
    
    
    #Step 2. populate nucl_elements with entity related elements - sides and histones (DNA?)
    nucl_elements={'side1':'(not all)','side2':'(not all)','H3':'(not all)','H4':'(not all)','H2A':'(not all)','H2B':'(not all)'}
    for k,v in components.items():
        if v['entity']=='histone' and v['type']!='H1':
            nucl_elements['side%d'%v['side']]=nucl_elements['side%d'%v['side']]+'or (segid %s)'%k
            nucl_elements[v['type']]=nucl_elements[v['type']]+'or (segid %s)'%k
    
    
    #Step 3. Populate with ss_elements
    for i in basic_elem_terms:
        nucl_elements[i]='(not all)'
        
    for k,v in components.items():
        if v['entity']=='histone':
            for k2,v2 in ss_seq_dict[k].items():
                nucl_elements[k2]=nucl_elements[k2]+' or (segid %s and resid %d:%d)'%(k,v2[0],v2[1])

    #Step 4. Populate with groups (super selections)
    for k,v in elem_groups.items():
        nucl_elements[k]='(not all)'
        for i in v:
            nucl_elements[k]=nucl_elements[k]+' or (%s)'%nucl_elements[i]
                
    
    #Step 5. Add some atom definitions of DNA
    nucl_elements['bases']='nucleic and (name N1 N2 N3 N4 N6 N7 N9 C2 C4 C5 C6 C7 C8 C4 O2 O6 C5M)'
    nucl_elements['sugars']='nucleic and (name O4\' C1\' C2\' C3\' C4\' C5\')'
    nucl_elements['phosphates']='nucleic and (name O5\' O3\' P OP1 OP2 O1P O2P)'
    nucl_elements['min_groove']='nucleic and ((resname DC and name C2 N3 O2) or (resname DG and name N2 C2 N3 C4) or (resname DA and name  C2 N3 C4) or (resname DT and name O2 C2 N3))'
    nucl_elements['cont']='nucleic and ((resname DC and name C6 C5 C4 N4) or (resname DG and name O6 C6 C5 N7 C8) or (resname DA and name  N6 C5 N7) or (resname DT and name C6 C5 C7 C5M O4 C4))'


    for k,v in components.items():
        if v['entity']=='histone' and v['type']!='H1':
            if 'nlp_id' in v.keys():
                npl_id='nlp_%s'%v['nlp_id']
                if not npl_id in nucl_elements.keys():
                    nucl_elements[npl_id]='(segid %s)'%k
                else:
                    nucl_elements[npl_id]=nucl_elements[npl_id][:-1]+' %s)'%k

    
    #as a precausion wrap everything in parenthesis
    for k,v in nucl_elements.items():
        nucl_elements[k]="(%s)"%v
        

    
    
    return nucl_elements


#old implementation
# def nucl_sel_expand(selection,elem_dict):
#     """
#     Given an element dict={'def':'resid 1:2'} and selection='name CA and def'
#     return 'name CA and resid 1:2
#     """
#     sel=selection
#     sel=' '+sel+' '
#     for k, v in elem_dict.items():
#         sel = sel.replace(' '+k+' ', ' '+v+' ')
#         sel = sel.replace('('+k, ' '+v+' ')
#         sel = sel.replace(k+')', ' '+v+' ')


    
#     return sel
#new implementation splist sel string with re and finds keywors one by one, may be a bit slower, but should be robust
import re
def nucl_sel_expand(selection,elem_dict):
    """
    Given an element dict={'def':'resid 1:2'} and selection='name CA and def'
    return 'name CA and resid 1:2
    """
    sel=selection
    # re search for word: [\w']+ or ( or ) or : or \s or ,  or . and create a list of search
    splitted_list=re.findall(r"[\w']+|\(|\)|:|\s+|,|.", sel)
    for index,item in enumerate(splitted_list):
        if item in elem_dict.keys():
            splitted_list[index]=elem_dict[item]
    return ''.join(splitted_list)
