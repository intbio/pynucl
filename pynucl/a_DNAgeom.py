"""
Analysis of DNA geometric parameters in NRF
relative twist
DNA path? - might be in a_geom also
deformation by the remodelers(?)


rTwist: we have C1' C1' vector (VEC, called BPOA in NAR paper) from top to bottom strand.
As a measure of rTw we can
v1)  measure angle between this vector and OXY plane.

v2) As in NAR Shaytan et al 2017 - take the angle between VEC and ORxy vector  (OR- is the vector from O to the mid-basepair, ORxy is its projection onto OXY plane, so its like a radial vector in cylindrical coordinate system(!))
this yielded almost identical results in NAR (see NAR where both v1 and v2 where compared). However, this has an advantage of being easier to calculate.

v3) Project the VEC to OZ-ORxy plane and calculate the angle between its projection and OXY plane.
this removes any components of VEC being deviating from OZ-OR plane, which it should follow in ideal superhelix.

NB OZ-ORxy plane is the same as OZ-OR plane

In this code v3 is implemented currently!!!!

An important question is the sign of the angle(!).


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
        self.bp_dict=nuclstr_instance.bp_dict
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
        
        
        # reimplemented using bp_dict
        top_resids=self.bp_dict['top'].keys()
        bot_resids=[self.bp_dict['top'][top_resid] for top_resid in top_resids]
        #lest shorten the selection too speed things up
        presel_top=self.u.select_atoms('nucleic and segid %s and name C1\''%nuclstr_instance.dyad_top[0])
        presel_bot=self.u.select_atoms('nucleic and segid %s and name C1\''%nuclstr_instance.dyad_bot[0])
        # be shure to leave all dublicates if altlocs
        presel_top=presel_top[(presel_top.altLocs=='') | (presel_top.altLocs=='A')]
        presel_bot=presel_bot[(presel_bot.altLocs=='') | (presel_bot.altLocs=='A')]
        
        # here we use mdanalysis ordered selections feature
        ordered_sel_list_top=['resnum %d'%resid for resid in top_resids]
        ordered_sel_list_bot=['resnum %d'%resid for resid in bot_resids]
        self.sel_top=presel_top.select_atoms(*ordered_sel_list_top)
        self.sel_bot=presel_bot.select_atoms(*ordered_sel_list_bot)
        
       
        
        # OLD implementation
        #self.sel_top=self.u.select_atoms('nucleic and segid %s and name C1\''%nuclstr_instance.dyad_top[0])
        #self.sel_bot=self.u.select_atoms('nucleic and segid %s and name C1\''%nuclstr_instance.dyad_bot[0])


        self.df_series=pd.DataFrame()

        for ts in tqdm(self.u.trajectory[self.time]):
            pos1=self.sel_top.positions
            pos2=self.sel_bot.positions # they are oredered with selection itself
            # OLD implementation
            #pos2=self.sel_bot.positions[::-1]
            vec=pos2-pos1 # THIS IS THE VEC
#             print(vec)
#             oz_vec=np.array([0,0,1])
            or_vec=(pos1+pos2)/2 # THIS IS OR vectore from center of nucleosme to center of base pair
            or_vec_xy=or_vec
            or_vec_xy[:,2]=0 # Here we project it to OXY plane
            or_vec_xy=or_vec_xy/norm(or_vec_xy,axis=1,keepdims=True) # And normalize
            
            #Now we will project vec onto OZ-ORxy plane and calulate angle there in the coordinate system formed by OZ and ORxy
            
            y=vec[:,2] # it's component along OZ will be now y
            
            x=np.sum(vec*or_vec_xy, axis=1) #it's x-component is length of its projection onto ORxy

            #alternatively as in NAR ??????? CHECK IT!!!
            #x=np.sign(np.dot(vec, or_vec_xy))*norm(np.array([vec[0],vec[1],0]))
            ang=np.arctan2(y, x)
            rTwf=np.degrees(ang) #-180 180
            rTw=np.abs(rTwf)
#             angle=np.degrees(np.arctan2(y, x))
#             rTwcont=np.degrees(np.unwrap(ang)) # this gives 2*pi noise if fisrt bb fluctuates a lot
            #this is a tricky way to unwrap starting from the dayd
            rTwcont=np.degrees(np.append(np.unwrap(ang[nuclstr_instance.DNA_left_length::-1])[-1:0:-1],np.unwrap(ang[nuclstr_instance.DNA_left_length:])))
            df=pd.DataFrame({'BPnum':range(1,len(self.bp_dict['top'].keys())+1),
                             'BPnum_dyad':np.array(list(self.bp_dict['top'].keys()))-nuclstr_instance.dyad_top[1],
                             'rTw':rTw,'rTwf':rTwf,'rTwcont':rTwcont})
            # OLD IMPL
            #df=pd.DataFrame({'BPnum':range(1,nuclstr_instance.DNA_top_length+1),'rTw':rTw,'rTwf':rTwf,'rTwcont':rTwcont})
            df['Frame']=ts.frame
            self.df_series=pd.concat([self.df_series,df])
        
        # OLD IMPL
        #self.df_series['BPnum_dyad']=self.df_series['BPnum']-nuclstr_instance.DNA_left_length-1
        self.df_series['segid']=nuclstr_instance.dyad_top[0]
        self.df_series['resid']=self.df_series['BPnum']-1+nuclstr_instance.seqs[nuclstr_instance.dyad_top[0]]['resid_start']
        
        self.df=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid'])['rTwcont'].apply(np.mean).reset_index()
        self.df['rTwf']= (self.df['rTwcont'] + 180) % (2 * 180) - 180
        self.df['rTw']= np.abs(self.df['rTwf'])
        self.df_std=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid'])['rTwcont'].agg('std').reset_index() 
        self.df_std['rTw']=self.df_std['rTwcont']
        self.df_std['rTwf']=self.df_std['rTwcont']




