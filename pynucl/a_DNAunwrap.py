"""
A hihgh level class to analyze DNA unwrapping
Analysis through contacts,
DNA deviation.

Plotting of the graphs.

"""

# import MDAnalysis as mda
# from MDAnalysis.analysis import align
# from MDAnalysis.analysis.rms import rmsd
from numpy.linalg import norm
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from .a_contacts import a_contacts
from .a_DNAparams import a_DNA
from .plot1d import plot_line_mpl
from .pynucl import nuclstr
import matplotlib.pyplot as plt

class a_DNAunw_cont:
    """
    Estimate DNA unwrapping from contact analsyis.
    The reference structure is the frame = 0.
    """
    
    def __init__(self,nuclstr_instance,time=None):
        
        #take
        #nuclstr_instance.u
        p=nuclstr_instance
        self.p=p
        if time is None:
            c=a_contacts(p,'nucleic','inner_core')
        else:
            c=a_contacts(p,'nucleic','inner_core',time=time)


        ds=c.get_num_int_profile_series()
        ds.loc[ds['segid']==p.DNA_bot_segid,'resid']=p.dyad_top[1]+p.dyad_bot[1]-ds.loc[ds['segid']==p.DNA_bot_segid,'resid']
        
        nds=ds.groupby(['resid','Frame']).agg('sum').reset_index()
        self.unw=nds.groupby(['Frame'])['resid'].agg(['min','max']).reset_index()
        self.unw['prox']=self.unw['min']-(p.dyad_top[1]-p.DNA_left_length)
        self.unw['dist']=p.dyad_top[1]+p.DNA_right_length-self.unw['max']





    def plot(self,ax=None,figsize=(10,4),dpi=300):
        if ax is None:
            fig,ax=plt.subplots(figsize=figsize,dpi=dpi)
        ax.plot(self.unw['Frame'],self.unw['prox'],label='Proximal end')
        ax.plot(self.unw['Frame'],self.unw['dist'],label='Distal end')
        ax.legend()
        ax.set_title('%s: Unwrapped base pairs, contact analysis'%self.p.name)
        ax.set_ylabel('Unwrapped base pairs')
        ax.set_xlabel('Time')
        return ax
    
    def plot_pos(self,ax=None,figsize=(10,4),dpi=300):
        if ax is None:
            fig,ax=plt.subplots(figsize=figsize,dpi=dpi)
        ax.plot(self.unw['Frame'],self.unw['min'],label='Proximal end')
        ax.plot(self.unw['Frame'],self.unw['max'],label='Distal end')
        ax.legend()
        ax.set_title('%s: Last wrapped DNA base pair, contact analysis'%self.p.name)
        ax.set_ylabel('DNA position')
        ax.set_xlabel('Time')
        return ax


class a_DNAunw_pos:
    """
    Estimate DNA unwrapping from position of DNA bp centers realative to reference
    The reference is the frame number or a structure (nuclstr_instance).
    """
    
    def __init__(self,nuclstr_instance,DNAparam=None,threshold=10,ref=0):
        

        p=nuclstr_instance
        self.p=p
        self.threshold=threshold
        if DNAparam is None:
            DNAparam=a_DNA(p)
        if isinstance(ref,nuclstr):
            ref=a_DNA(ref).df
        DNAparam.calc_DNA_bp_closest(ref=ref)
        df=DNAparam.df_series.copy()
        ds=df.loc[df['BP_mindist']<threshold]
        
        ds.loc[ds['segid']==p.DNA_bot_segid,'resid']=p.dyad_top[1]+p.dyad_bot[1]-ds.loc[ds['segid']==p.DNA_bot_segid,'resid']
        nds=ds.groupby(['resid','Frame']).agg('sum').reset_index()
        self.unw=nds.groupby(['Frame'])['resid'].agg(['min','max']).reset_index()
        self.unw['prox']=self.unw['min']-(p.dyad_top[1]-p.DNA_left_length)
        self.unw['dist']=p.dyad_top[1]+p.DNA_right_length-self.unw['max']
        

    def plot_pos(self,ax=None,figsize=(10,4),dpi=300):
        if ax is None:
            fig,ax=plt.subplots(figsize=figsize,dpi=dpi)
        ax.plot(self.unw['Frame'],self.unw['min'],label='Proximal end')
        ax.plot(self.unw['Frame'],self.unw['max'],label='Distal end')
        ax.legend()
        ax.set_title('%s: Last wrapped DNA base pair, distance analysis analysis'%self.p.name)
        ax.set_ylabel('DNA position')
        ax.set_xlabel('Time')
        return ax

    def plot(self,ax=None,figsize=(10,4),dpi=300):
        if ax is None:
            fig,ax=plt.subplots(figsize=figsize,dpi=dpi)
        ax.plot(self.unw['Frame'],self.unw['prox'],label='Proximal end')
        ax.plot(self.unw['Frame'],self.unw['dist'],label='Distal end')
        ax.legend()
        ax.set_title('%s: Unwrapped base pairs, distance analysis'%self.p.name)
        ax.set_ylabel('Unwrapped base pairs')
        ax.set_xlabel('Time')
        return ax
    
        
    def plot_fluct(self,ax=None,figsize=(10,5),dpi=300):
        if ax is None:
            fig,ax=plt.subplots(figsize=figsize,dpi=dpi)
            ax=plot_line_mpl(DNAparam.df_series,self.p,column='BP_mindist',\
                                    ymin=0,ref=None,startnumber=p.dyad_top-p.DNA_left_length,\
                             funcgroups='\\funcgroup{xxx}{A}{Black}{Green}{upper}{up}',color_ref='black',color='#AAAAAA')
            dump=ax.set_xticks([4,14,24,34,44,54,64,74,84,94,104,114,124,134])
            dump=ax.set_title('%s: Deviation of DNA from reference path'%self.p.name)
            dump=ax.set_ylabel('Max.dev.,A.')
            ax.set_axisbelow(False)
            ax.grid(axis='x', color='0.0',zorder=3)
            ax.set_xlabel('DNA position')
        return ax