# -*- coding: utf-8 -*-
#This module provides routines to identify positions of histone features/elements (helices, loops, sheets, etc.) within a histone sequence.

#Specifically:
# 1) This module provides a vocabulary for template sequences (1KX5 and a few others)
# 2) Routines to transfer this vocubulary to other histones via alignment.
# 3) Generation of shading features for pytexshade is also supported.


# The elements are as in latest update of HistoneDB 
# https://github.com/intbio/histonedb/blob/master/static/browse/info/features.json
# see also https://github.com/intbio/seq_tools/blob/master/seq_tools/hist_ss.py (details of docking domain might differ there)


###Element definition protocol
# 1) Main reference structure is 1KX5 with a look at 1AOI
# 2) Helices are defined by taking the minimum length of alpha helices over symmetric chains in 1kx5
# 3) docking domain thoughts:
# The  differences may rise in definition of H2A-docking domain:
# In Luger 1997 1aoi paper  it is defined as 80 – 119 (1 based numbering)
# In Suto Luger 2000 H2A.Z paper  81-119
# Shaytan et al. JMB 2016 Fig.2. 80-118
# Original HistoneDB 81-119
# In the latest update of HistoneDB we thought to use Luger 1997 convention 80-119, but now we will use 80-118
# However, in MD we see that K119 is too flexible, so it makes sense to define docking domain as 80-118,
# as we do here

#H2A.Z mapping should be provided separately, manual structural alignment used to precisely place the loop elements.


# Numbering 
# We use the sequences that correspond to those present in the PDB.
# in case of 1kx5 these are identical to SEQREC, except for H2B where SEQREC has a 3 aa ovehang at N-terminus
# On those sequences we map the elements via 0-based numbering, i.e. they are one off from resid.
# 'element_name':[start,stop], start stop - are included in the range (this is not as in Python), this is due to historical compatibility with TexShade and resid start:stop selection in MDanalysis.


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import pairwise2
from Bio import Align
from Bio.Align import substitution_matrices

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

#1KX5 seq for reference, these are sequences of resolved chains not SEQREC (with seqrec h2b there is a problem it has a 3 aa overhang on N-tail)

templ_H3 =Seq("ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVALFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA", IUPAC.protein)
templ_H4 = Seq("SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG", IUPAC.protein)
templ_H2A = Seq("SGRGKQGGKTRAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLPKKTESSKSKSK", IUPAC.protein)
templ_H2AZ = Seq("AGGKAGKDSGKAKTKAVSRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIGKKGQQKTV", IUPAC.protein) #from 1f66
templ_H2B = Seq("AKSAPAPKKGSKKAVTKTQKKDGKKRRKTRKESYAIYVYKVLKQVHPDTGISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK", IUPAC.protein)


#'element_name':[start,stop], start stop - are both in the range


#core - are parts that we modelled as 1kx5_ntm
#inner_core - tails truncated further to alphaN or alpha1 tail boders.

ss_templ_H3={'inner_core':[43,134],'core':[38,134],'alphaN':[43,56],'alpha1':[62,76],'alpha2':[84,113],'alpha3':[119,130],'loopL1':[78,83],'loopL2':[114,118],'beta1':[82,83],'beta2':[117,118],'mgarg1':[62,62],'mgarg2':[82,82],'mgarg3':[48,48]}
ss_templ_H4={'inner_core':[23,101],'core':[24,101],'alpha1ext':[23,28],'alpha1':[29,40],'alpha2':[48,75],'alpha3':[81,92],'loopL1':[41,47],'loopL2':[76,81],'beta1':[44,45],'beta2':[79,80],'beta3':[95,97],'mgarg1':[44,44]}
ss_templ_H2A={'inner_core':[15,117],'core':[15,117],'alpha1ext':[15,21],'alpha1':[25,36],'alpha2':[45,72],'alpha3':[78,88],'alpha3ext':[89,96],'loopL1':[37,44],'loopL2':[73,77],'beta1':[41,42],'beta2':[76,77],'beta3':[99,101],'docking_domain':[79,117],'mgarg1':[41,41],'mgarg2':[76,76]}
ss_templ_H2B={'inner_core':[33,121],'core':[25,121],'alpha1':[33,45],'alpha2':[51,80],'alpha3':[86,98],'alphaC':[99,119],'loopL1':[46,50],'loopL2':[81,85],'beta1':[49,50],'beta2':[84,85],'mgarg1':[29,29]}

ss_templ_H2AZ={'alpha1ext':[17,23],'alpha1':[27,38],'alpha2':[48,75],'alpha3':[81,91],'alpha3ext':[92,99],'loopL1':[39,47],'loopL2':[76,80],'beta1':[44,45],'beta2':[79,80],'beta3':[101,103],'docking_domain':[82,119],'mgarg1':[44,44],'mgarg2':[79,79]}

ftypes_dict={'alpha1ext':'helix','alphaN':'helix','alpha1':'helix','alpha2':'helix','alpha3':'helix','beta1':'sheet','beta2':'sheet','beta3':'sheet','loopL2':'loop','loopL1':'loop','mgarg1':'frameblock','mgarg2':'frameblock','mgarg3':'frameblock'}
ftext_dict={'alpha1':'$\\alpha 1$','alpha2':'$\\alpha 2$','alpha3':'$\\alpha 3$','alphaN':'$\\alpha N$','beta1':'$\\beta 1$','beta2':'$\\beta 2$'}

ss_templ_dict={'H4':{'seq':templ_H4,'elem':ss_templ_H4},
            'H3':{'seq':templ_H3,'elem':ss_templ_H3},
            'H2A':{'seq':templ_H2A,'elem':ss_templ_H2A},
            'H2A.Z':{'seq':templ_H2AZ,'elem':ss_templ_H2AZ},
            'H2B':{'seq':templ_H2B,'elem':ss_templ_H2B}}


basic_elem_terms=set(ss_templ_H3.keys()).union(set(ss_templ_H4.keys())).union(set(ss_templ_H2A.keys())).union(set(ss_templ_H2B.keys()))

elem_groups={'hist_folds':['alpha1','alpha2','alpha3'],
                 'mgargs':['mgarg1','mgarg2','mgarg3'],
                 'loops':['loopL1','loopL2']}


def hist_ident(seq,hist_type=None):
    """
    identifies the best matching sequence from ss_templ_dict
    """
    best_score=100
    best_seq='undef'
    best_aln='no'
    if hist_type is None:
        d=ss_templ_dict
    else:
        d={'hist_type':ss_templ_dict['hist_type']}
    for k,v in d.items():
        i=0
        for a in aligner.align(str(v['seq']),str(seq)):
            if(a.score>best_score):
                best_score=a.score
                best_seq=k
                best_aln=a
            i=i+1
            if i>100:
                break
#     print(best_score)
    return best_seq,best_aln
            

    
def maploc(aln,tstart,tstop):
    t2q={}
    
    for i,j in zip(aln.aligned[0],aln.aligned[1]):
        for x,y in zip(range(*i),range(*j)):
            t2q[x]=y
            
    for i in range(len(aln.target)):
        qstart=t2q.get(tstart,-1)
        if qstart>=0:
            break
        else:
            tstart=tstart+1

    for i in range(len(aln.target)):
        qstop=t2q.get(tstop,-1)
        if qstop>=0:
            break
        else:
            tstop=tstop-1

    if qstart<=qstop:
        return qstart,qstop
    else:
        return -1,-1
    

def hist_features(seq,hist_type=None):
    """
    should return a list of secondary structure elements like ss_templ_* dicts but for a specific sequence, by aligning it to the template (specified by type or idendified by ident_hist)
    """
    hist_type,aln=hist_ident(seq)
    if hist_type=='undef':
        return []
    feat=[]
    for k,v in ss_templ_dict[hist_type]['elem'].items():
        start,stop=maploc(aln,v[0],v[1])
        if start>=0 and stop>=0:
            feat.append(SeqFeature(FeatureLocation(start,stop+1), type=ftypes_dict.get(k,'misc'),id=k))
    return feat
    

def hist_shade_features(input_features,force_feature_pos='bottom'):
    """
    Generate a fetaures list for pytexshade for histones from list of biopython SeqFeature objects
    """
    fc={'sheet':'-->'}
    #In texshade [1,1] will be colored as one residue. this is different from biopython where [1,1] selects nothing!
    #create a feature overlap string
    features=[]
    for f in input_features:
        pos=force_feature_pos if force_feature_pos else 'top'
        if(f.id=='core'):
            pos=pos[0]+pos[0]+pos
        if(f.type=='loop'):
            pos=pos[0]+pos
        features.append({'style':fc.get(f.type,f.type),'seqref':'1','sel':[int(f.location.start),f.location.end-1],'position':pos,'text':ftext_dict.get(f.id,f.id)})

    return features

def hist_shf4seq(seq):
    """
    for a given histone sequence will return a list of features for shading with pytexshade
    """
    return hist_shade_features(hist_features(seq))

    
            # Our pytexshade utils take annotation in their own format
        # However their is a converter from Biopython-style features in SecRec class to pytexshade (shade.seqfeat3shadefeat), and an example of getting annotations from NCBI (seqplot/pdb_plot.py)
        # but the type of features is limited there. So for histones we stick with our format.
        # we will need following
        # features=[{'style':'-->','helix','loop','frameblock'},'sel':[begin(0-based),end(inclusive)],'text':'\alpha 1']
