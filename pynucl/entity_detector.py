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
from .seq_utils import seqrec_from_pdb,seq_from_mda

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
            logger.debug('No top and bot strands detected from database, finding two complementary chains and assigning by chain name')
            dna_chain_combinations=list(itertools.combinations(dna_chain_list,2))
            score=[]
            for segid1,segid2 in dna_chain_combinations:
                score.append(get_chains_complementarity(entity_dict[segid1]['sequence'],entity_dict[segid2]['sequence']))
            if max(score)>0.85:
                top_segid,bot_segid=dna_chain_combinations[score.index(max(score))]
                logger.debug(f'Chains {top_segid} and {bot_segid} are complementary assigning top and bot by chain name')
        if not(top_segid is None) and not (bot_segid is None):
            entity_dict[top_segid]['type']='DNAtop'
            entity_dict[bot_segid]['type']='DNAbot'
        else:
            logger.debug('Could not find chains with significant complementarity')
    # TODO check that the sewuence detected is the same )
    return (entity_dict)
    
    
    
    
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




