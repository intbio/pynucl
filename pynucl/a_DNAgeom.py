"""
Analysis of DNA geometric parameters in NRF
relative twist
DNA path? - might be in a_geom also
deformation by the remodelers(?)


rTwist: we have C1' C1' vector (VEC) from top to bottom strand.
As a measure of rTw:
we would like to measure angle between OXY plane.
or alternatively we can first project the VEC
is the OZ OR angle, where OZ - is the z-axis, OR - is the axis from O to the mid-basepair

The OZ-OR angle is the angle between VEC and the plane of OZOR 

"""

# import MDAnalysis as mda
# from MDAnalysis.analysis import align
# from MDAnalysis.analysis.rms import rmsd
from numpy.linalg import norm
import numpy as np
import pandas as pd
from tqdm.auto import tqdm


class a_DNArtw:
    """
    Get DNA relative twist values
    """
    
    def __init__(self,nuclstr_instance,time=None):
        
        #take
        #nuclstr_instance.u
        self.u=nuclstr_instance.u
        if time is None:
            if isinstance(nuclstr_instance.time,slice):
                self.time=nuclstr_instance.time
            else:
                self.time=slice(*nuclstr_instance.time)
        else:
            if isinstance(time,slice):
                self.time=time
            else:
                self.time=slice(*time)
           

        self.sel_top=self.u.select_atoms('nucleic and segid %s and name C1\''%nuclstr_instance.dyad_top[0])
        self.sel_bot=self.u.select_atoms('nucleic and segid %s and name C1\''%nuclstr_instance.dyad_bot[0])


        self.df_series=pd.DataFrame()

        for ts in tqdm(self.u.trajectory[self.time]):
            pos1=self.sel_top.positions
            pos2=self.sel_bot.positions[::-1]
            vec=pos2-pos1
#             print(vec)
#             oz_vec=np.array([0,0,1])
            or_vec=(pos1+pos2)/2
            or_vec[:,2]=0
            or_vec=or_vec/norm(or_vec,axis=1,keepdims=True)
            y=vec[:,2]
            x=np.sum(vec*or_vec, axis=1)
#             x = np.dot(vec, or_vec)

            #alternatively as in NAR
            #x=np.sign(np.dot(vec, or_vec))*norm(np.array([vec[0],vec[1],0]))
            
            rTwf=np.degrees(np.arctan2(y, x))
            rTw=np.abs(rTwf)
#             angle=np.degrees(np.arctan2(y, x))
            
            df=pd.DataFrame({'BPnum':range(1,nuclstr_instance.DNA_top_length+1),'rTw':rTw,'rTwf':rTwf})
            df['Time']=ts.frame
            self.df_series=pd.concat([self.df_series,df])
        
        
        self.df_series['BPnum_dyad']=self.df_series['BPnum']-nuclstr_instance.DNA_left_length-1
        self.df_series['segid']=nuclstr_instance.dyad_top[0]
        self.df_series['resid']=self.df_series['BPnum']-1+nuclstr_instance.seqs[nuclstr_instance.dyad_top[0]]['resid_start']
        
        self.df=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid'])['rTw','rTwf'].apply(np.mean).reset_index()
        self.df_std=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid'])['rTw','rTwf'].agg('std').reset_index() 


