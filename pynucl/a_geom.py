"""
Analysis of system's geometry

# Scalar values
1) Distances between centers-of-mass of two selections (as a reduced case - between two the atoms)
2) Angles for 3 selections
3) Dihedrals for 4 selections
4) RMSD of a group with respect to ref.
5) Angle between helices: 
   reference helices are aligned and vectors defined as
   connecting first and last CA atoms

# Vector values
1) Coordinates of centers-of-mass of a selection
2) Vectors between two selections
3) Vectors of a helix, via first fitting a reference helix

# Multiplexed input/output for a list of selections, or demultiplexed by residue

# Matrix data
1)?correlations for PCA?

# The logic of combining one structure and trj-analysis:
always return average, and have also a _series object. (see pymolint for approx. reference).
#Also we will have rmsf functions.

# The expected input:
sel1 or (sel1, sel2, ...) - tuple strictly, NOT LIST. Depends if we need one or many selections to calculate the geometric parameter.
## For multiplexed input
[sel1,sel2] or [(sel1,sel2, ...),(sel1,sel2, ... )]
## Alternatively, specify split_by_resid=True
This will demultiplex spections by segid,resid

# The expected output format is:
## A dataframe for average
value
 0.5
## A dataframe for series
Time, value
1, 0.5
# For multiplexed input/output
  value
1  0.5
2  0.6

"""

# from pymolint.mol_int import struct2cont # this is for contact analysis
from .nucl_meta_select import create_elem_dict,nucl_sel_expand
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from numpy.linalg import norm
from tqdm.auto import tqdm




class mol_geom:
    """
    Class to analyze geometry in  molecular systems
    See ideas at the top of the file.
    types of mol_geometry:
       --scalar--
        dist : selections= ('resid 5','resid 6') or [(,),(,)]
        ang  : selections= ('resid 5','resid 6','resid 7') or [(,,),(,,)..]      
        dih  : selections= ('','','','') or [(,,,),(,,,,)...] 
        rmsd : selections= 'resid 5' or [(),()]
        helixang : selections= ('helix1','helix2') or [(,),(,)]
        
       --vector--
       coord : selections= 'resid 5' or [(),()]
       vec   : selections= ('resid 5','resid 6','resid 7') or [(,,),(,,)..] 
       helixvec : selections= 'helix1' or [(),()]
       
       Inputs:
       sel_split=False or str - demultiplex by 'resid' or 'atom'
          only for coord, rmsd, dist
       ref  - reference for rmsd or helix* parameters calculations.
              May be a Universe.
              or a frame number,
              by default is first frame (0 frame).
       rmsd_superpostion - superposition parameter to rmsd.
      
    """
    def __init__(self,u,geom_type,selections,ref=0,sel_split=False, time=(None,None,None), rmsd_superposition=False,resnumfix=False, **kwargs):
        self.u=u
        self.geom_type=geom_type
        self.ref=ref
        self.sel_split=sel_split
        self.rmsd_superposition=rmsd_superposition
        self.resnumfix=resnumfix
        
        #Let bring selections to a unified system:
        #Always a list of tuples.
        if isinstance(selections,list):
            self.selections=[]
            for i in selections:
                if isinstance(i,tuple):
                    self.selections.append(i)
                elif isinstance(i,str):
                    self.selections.append((i,))
                else:
                    raise Exception('Bad selection format')
        if isinstance(selections,tuple):
            self.selections=[selections]
        if isinstance(selections,str):
            self.selections=[(selections,)]
            
        
#         #########make real selections##########
        self.s=self.sels_text2obj(self.u)
        if isinstance(self.ref,mda.Universe):
            self.refs=self.sels_text2obj(self.ref)
            
        self.time=list(time)
        
        if time[0] is None:
            self.time[0]=0
        if time[1] is None:
            self.time[1]=self.u.trajectory.n_frames
        if time[2] is None:
            self.time[2]=1
        
        #Scalars
        if self.geom_type=='dist':
            self.dist()
        if self.geom_type=='ang':
            self.ang()
        if self.geom_type=='dih':
            self.dih()
        if self.geom_type=='rmsd':
            self.rmsd()
#         if self.type=='ang_btw_hel': #TODO - see mol_geom
#             self.ang_btw_hel()


        #Vectors
        if self.geom_type=='coord':
            self.coord()



#         if self.geom_type=='helix_vector': #TODO - see mol_geom
#             self.helix_vector()

##### Block to get annotated dataframes and get some stats

        self.df_series=pd.DataFrame({'Frame':[i for i in range(self.time[0],self.time[1],self.time[2]) for j in range(len(self.s))],'sel_num':[j for i in range(self.time[0],self.time[1],self.time[2]) for j in range(len(self.s))],self.geom_type:self.values})
                
        self.df=self.df_series.groupby(['sel_num'])[self.geom_type].apply(np.mean).reset_index()
        self.df_std=self.df_series.groupby(['sel_num'])[self.geom_type].apply(lambda x: np.std(x.to_numpy(),ddof=1)).reset_index() 
        
        
        if self.geom_type=='coord':
            #we will also expand to x,y,z for compatibility
            self.df_series['x']=[i[0] for i in self.df_series['coord']]
            self.df_series['y']=[i[1] for i in self.df_series['coord']]
            self.df_series['z']=[i[2] for i in self.df_series['coord']]
            self.df['x']=[i[0] for i in self.df['coord']]
            self.df['y']=[i[1] for i in self.df['coord']]
            self.df['z']=[i[2] for i in self.df['coord']]
            self.df_std['x']=[i[0] for i in self.df_std['coord']]
            self.df_std['y']=[i[1] for i in self.df_std['coord']]
            self.df_std['z']=[i[2] for i in self.df_std['coord']]
        
        
        
        if sel_split=='resid' or sel_split=='atom':#add info to the dataframe
            self.df_series['segid']=np.tile(self.ds.segids,len(range(self.time[0],self.time[1],self.time[2])))
            self.df['segid']=self.ds.segids
            self.df_std['segid']=self.ds.segids

            self.df_series['resid']=np.tile(self.ds.resids,len(range(self.time[0],self.time[1],self.time[2])))
            self.df['resid']=self.ds.resids
            self.df_std['resid']=self.ds.resids

        if sel_split=='atom':#add info to the dataframe
            self.df_series['name']=np.tile(self.ds.names,len(range(self.time[0],self.time[1],self.time[2])))
            self.df['name']=self.ds.names
            self.df_std['name']=self.ds.names

            self.df_series['index']=np.tile(self.ds.indices,len(range(self.time[0],self.time[1],self.time[2])))
            self.df['index']=self.ds.indices
            self.df_std['index']=self.ds.indices    

#######
    
    def sels_text2obj(self,u):
        selobj=[]
        if self.sel_split==False:
            for i in self.selections:
                selobj.append([u.select_atoms(j) for j in i])
        elif self.sel_split=='resid': #demultiplex first sel[0][0] in a list of selections by resid
            self.ds=u.select_atoms(self.selections[0][0]).residues
            for si,ri in zip(self.ds.segids,self.ds.resids):
                if self.resnumfix:
                    selobj.append([u.select_atoms('segid %s and resnum %s and (%s)'%(si,ri,self.selections[0][0]))])
                else:
                    selobj.append([u.select_atoms('segid %s and resid %s and (%s)'%(si,ri,self.selections[0][0]))])


        elif self.sel_split=='atom': #demultiplex first sel[0][0] in a list of selections by atoms
            self.ds=u.select_atoms(self.selections[0][0]).atoms
            for ai in self.ds.indices:
                selobj.append([u.select_atoms('index %d'%(ai))])
        return selobj
    
    
    
    def dist(self):
        self.values=[] #first cycle is  over time, second  is multipex
        for ts in tqdm(self.u.trajectory[self.time[0]:self.time[1]:self.time[2]]):
            self.values.extend([norm(i[0].center_of_mass() - i[1].center_of_mass()) for i in self.s])    
     
    def rmsd(self):
        self.values=[] #first cycle is  over time, second  is multipex
        if isinstance(self.ref,int):
            self.u.trajectory[self.ref]
            refpos=[i[0].positions for i in self.s]
            self.u.trajectory[0]
        if isinstance(self.ref,mda.Universe):
            refpos=[i[0].positions for i in self.refs]

        for ts in tqdm(self.u.trajectory[self.time[0]:self.time[1]:self.time[2]]):
            self.values.extend([rmsd(r,i[0].positions, superposition=self.rmsd_superposition) for r,i in zip(refpos,self.s)]) 
    
    
    def ang(self):
        self.values=[]
        for ts in tqdm(self.u.trajectory[self.time[0]:self.time[1]:self.time[2]]):
            snapvals=[]
            for i in self.s:
                v1=i[0].center_of_mass()-i[1].center_of_mass()
                v2=i[2].center_of_mass()-i[1].center_of_mass()
                cosang = np.dot(v1, v2)
                sinang = norm(np.cross(v1, v2))
                snapvals.append(np.degrees(np.arctan2(sinang, cosang)))
            self.values.extend(snapvals)
    
    #TODO that one shoul be near instant with u.trajectory.timeseries
    def coord(self):
        self.values=[] #first cycle is over time, second  is multipex
        for ts in tqdm(self.u.trajectory[self.time[0]:self.time[1]:self.time[2]]):
            self.values.extend([i[0].center_of_mass() for i in self.s]) 
        
    def dih(self):
        self.values=[] #first cycle is over time, second  is multipex
        for ts in tqdm(self.u.trajectory[self.time[0]:self.time[1]:self.time[2]]):
            snapvals=[]
            for i in self.s:
                snapvals.append(dihedral([i[0].center_of_mass(),i[1].center_of_mass(),i[2].center_of_mass(),i[3].center_of_mass()]))
            self.values.extend(snapvals)
        
        
        
class a_geom(mol_geom):
    """
    Class to analyze geometry parameters in nucleosome
    Get two selections (using expanded meta synataxis)
    """
    def __init__(self,nuclstr_instance,geom_type,selections, **kwargs):
        if 'time' in kwargs:
            pass
        else:
            kwargs['time']=nuclstr_instance.time
            
        #will need to expand selections for every possible input nested list/tuple input format of selections    
        def rec_sel_update(s,expdict):
            if isinstance(s,str):
                return nucl_sel_expand(s,expdict)
            if isinstance(s,tuple):
                return tuple([rec_sel_update(i,expdict) for i in s])
            if isinstance(s,list):
                return [rec_sel_update(i,expdict) for i in s]
            
        sels=rec_sel_update(selections,nuclstr_instance.nucl_elements)
        
        super().__init__(nuclstr_instance.u,geom_type,sels, **kwargs)


def dihedral(p):
    """
    Taken from here
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


