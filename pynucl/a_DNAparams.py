"""
Streamlining calculation of DNA parameters by 3DNA (CURVES in the future)
for pynucl objects 


TODO: we need to make it consederably faster - without dumping multitude of files
"""
from io import StringIO
from tqdm.auto import tqdm
from multiprocessing import Pool
import DNAtools
import tempfile
import os
import uuid
import pandas as pd
import numpy as np 
import logging
from MDAnalysis import Merge
import scipy
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)s %(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

class a_DNA:
    """
    Analyze DNA parameters for a pynucl object
    bps - define what basepairs to consider in calculating params.
        None - all
        [] - list of resids along the top chain
    Will provide df, df_series, df_std
    """
    
    def __init__(self,nuclstr_instance,time=None,bps=None,num_threads=1,dir=None):
        
        #take
        #nuclstr_instance.u
        self.resid_start=nuclstr_instance.seqs[nuclstr_instance.dyad_top[0]]['resid_start']
        self.DNA_length=nuclstr_instance.DNA_top_length
        self.bp_dict=nuclstr_instance.bp_dict
        if bps is None:
            #MAJOR CHANGES - idea is that we use bp_dict object, which is derived at structure initialisation
            # using allignment with rev compl bot strand
            #self.bps=np.arange(self.resid_start,self.resid_start+self.DNA_length)
            self.bps=np.array(list(nuclstr_instance.bp_dict['top'].keys()))
        else:
            logger.warning('With bps set manually, sugar information is incorrect!!!')
            self.bps=[]
            for bp in bps:
                if bp in (nuclstr_instance.bp_dict['top'].keys()):
                    self.bps.append(bp)
            self.bps=np.array(self.bps)
        self.p=nuclstr_instance
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
           
#         print(self.time)
        self.sel=self.u.select_atoms('nucleic') #should be allright now
        self.sel_u=Merge(self.sel)
        #OLD WAY
#         with tempfile.TemporaryDirectory() as TEMP:
#             logger.debug('created temporary directory ' + TEMP)
#             #At first we need to makup a couple of unique file names
#             self.u.trajectory[self.time.start]
#             self.sel.write(os.path.join(TEMP,'first.pdb'))
#             self.bp_ref_file=DNAtools.X3DNA_find_pair(os.path.join(TEMP,'first.pdb'))

        
        #####New way
        self.bp_ref_file=self.gen_find_pair_file()
        
        
        
        #Run X3DNA_analyze
        if num_threads!=1:
                            
            with tempfile.TemporaryDirectory(dir=dir) as TEMP:
#                 print(TEMP)
                pdbfiles=[]
                frames=[]
                for ts in self.u.trajectory[self.time]:
                    unique=str(uuid.uuid4())
                    pdb = unique+'.pdb'
                    pdbfiles.append(os.path.join(TEMP,pdb))
                    frames.append(ts.frame)
                with Pool(num_threads) as pool:
                    # to be picklable we must construct zip with such an order pdb,frame,sel,bp_ref_file,bps
                    N=len(pdbfiles)
                    df_list = list(tqdm(pool.imap(X3DNA_analyze_runner, zip(pdbfiles,
                                                                            frames,
                                                                            [self.sel]*N,
                                                                            [self.bp_ref_file]*N,
                                                                            [self.bps]*N)
                                                 ),
                                        total=len(pdbfiles)
                                       )
                                  )
                for df in df_list:
                    df['BPnum_dyad']=self.bps-self.p.dyad_top[1]
                    df['segid']=nuclstr_instance.dyad_top[0]
                    
                self.df_series=pd.concat(df_list)

#old singlethread
        elif num_threads==1:
            self.df_series=pd.DataFrame()
            with tempfile.TemporaryDirectory(dir=dir) as TEMP:
                for ts in tqdm(self.u.trajectory[self.time]):
                    unique=str(uuid.uuid4())
                    pdb = unique+'.pdb'
                    self.sel.write(os.path.join(TEMP,pdb))
                    df=DNAtools.X3DNA_analyze(os.path.join(TEMP,pdb),ref_fp_id=self.bp_ref_file)
    #                 print(df)
                    df=df.iloc[0:len(self.bps)]
                    df['Frame']=ts.frame
#                     df['Time']=ts.time
                    df['BPnum_dyad']=self.bps-self.p.dyad_top[1]
                    df['segid']=nuclstr_instance.dyad_top[0]
                    df['resid']=self.bps
                    self.df_series=pd.concat([self.df_series,df])
        
        
        #TODO: Add tricks if top and bottom in PDB are not in the the order top,bottom, but bottom,top

        self.df_series=self.df_series.reset_index(drop=True)
        
        self.df_series.replace('---',np.nan,inplace=True)
        self.df_series['alpha_1']=self.df_series['alpha_1'].astype(float)
        self.df_series['beta_1']=self.df_series['beta_1'].astype(float)
        self.df_series['epsilon_1']=self.df_series['epsilon_1'].astype(float)
        self.df_series['zeta_1']=self.df_series['zeta_1'].astype(float)
        self.df_series['epsilon_2']=self.df_series['epsilon_2'].astype(float)        
        self.df_series['e-z_1']=self.df_series['e-z_1'].astype(float)        
        self.df_series['alpha_2']=self.df_series['alpha_2'].astype(float)
        self.df_series['beta_2']=self.df_series['beta_2'].astype(float)
        self.df_series['zeta_1']=self.df_series['zeta_1'].astype(float)
        self.df_series['e-z_2']=self.df_series['e-z_2'].astype(float)  
        self.df_series['zeta_2']=self.df_series['zeta_2'].astype(float)        
        self.df_series['ssZp_1']=self.df_series['ssZp_1'].astype(float)
        self.df_series['ssZp_2']=self.df_series['ssZp_2'].astype(float)
        self.df_series['Dp_1']=self.df_series['Dp_1'].astype(float)
        self.df_series['Dp_2']=self.df_series['Dp_2'].astype(float)        
        
        
        try:
            self.do_stats()
        except:
            logger.warning('statistics failed!')
        
    def do_stats(self):
        # old version does not work
        #self.df=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid'])['x','y','z','Shear','Stretch','Stagger','Buckle','Prop-Tw','Opening','Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist','Pairing', 'chi_1', 'alpha_1', 'beta_1', 'gamma_1', 'delta_1','epsilon_1', 'zeta_1', 'e-z_1', 'chi_2', 'alpha_2','beta_2', 'gamma_2', 'delta_2', 'epsilon_2', 'zeta_2', 'e-z_2',  'v0_1', 'v1_1', 'v2_1', 'v3_1', 'v4_1', 'tm_1', 'P_1', 'ssZp_1', 'Dp_1', 'v0_2', 'v1_2', 'v2_2', 'v3_2', 'v4_2','tm_2', 'P_2', 'ssZp_2', 'Dp_2'].apply(np.mean).reset_index()
        #self.df=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid']).apply(np.mean).drop(['BPnum','BPnum_dyad','segid','resid'],level=4).unstack().reset_index()
        self.df=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid']).agg('mean').drop(columns=['Frame']).reset_index()
        
        # refactor
        #self.df_std=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid'])['x','y','z','Shear','Stretch','Stagger','Buckle','Prop-Tw','Opening','Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist','Pairing', 'chi_1', 'alpha_1', 'beta_1', 'gamma_1', 'delta_1','epsilon_1', 'zeta_1', 'e-z_1', 'chi_2', 'alpha_2','beta_2', 'gamma_2', 'delta_2', 'epsilon_2', 'zeta_2', 'e-z_2',  'v0_1', 'v1_1', 'v2_1', 'v3_1', 'v4_1', 'tm_1', 'P_1', 'ssZp_1', 'Dp_1', 'v0_2', 'v1_2', 'v2_2', 'v3_2', 'v4_2','tm_2', 'P_2', 'ssZp_2', 'Dp_2'].agg('std').reset_index()
        self.df_std=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid']).agg('std').drop(columns=['Frame']).reset_index()
        
        #self.df_min=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid'])['x','y','z','Shear','Stretch','Stagger','Buckle','Prop-Tw','Opening','Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist','Pairing', 'chi_1', 'alpha_1', 'beta_1', 'gamma_1', 'delta_1','epsilon_1', 'zeta_1', 'e-z_1', 'chi_2', 'alpha_2','beta_2', 'gamma_2', 'delta_2', 'epsilon_2', 'zeta_2', 'e-z_2',  'v0_1', 'v1_1', 'v2_1', 'v3_1', 'v4_1', 'tm_1', 'P_1', 'ssZp_1', 'Dp_1', 'v0_2', 'v1_2', 'v2_2', 'v3_2', 'v4_2','tm_2', 'P_2', 'ssZp_2', 'Dp_2'].agg('min').reset_index()
        self.df_min=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid']).agg('min').drop(columns=['Frame']).reset_index()
        
        #self.df_max=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid'])['x','y','z','Shear','Stretch','Stagger','Buckle','Prop-Tw','Opening','Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist','Pairing', 'chi_1', 'alpha_1', 'beta_1', 'gamma_1', 'delta_1','epsilon_1', 'zeta_1', 'e-z_1', 'chi_2', 'alpha_2','beta_2', 'gamma_2', 'delta_2', 'epsilon_2', 'zeta_2', 'e-z_2',  'v0_1', 'v1_1', 'v2_1', 'v3_1', 'v4_1', 'tm_1', 'P_1', 'ssZp_1', 'Dp_1', 'v0_2', 'v1_2', 'v2_2', 'v3_2', 'v4_2','tm_2', 'P_2', 'ssZp_2', 'Dp_2'].agg('max').reset_index()
        self.df_max=self.df_series.groupby(['BPnum','BPnum_dyad','segid','resid']).agg('max').drop(columns=['Frame']).reset_index()
       
        #for what ever reason .apply(lambda x: np.std(x.to_numpy(),ddof=1)).reset_index()  does not work(!)
        #Let's group x y z as coord array for compatibility
        self.df_series['coord']=[np.array([x,y,z]) for x,y,z in zip(self.df_series['x'],self.df_series['y'],self.df_series['z'])]
        self.df['coord']=[np.array([x,y,z]) for x,y,z in zip(self.df['x'],self.df['y'],self.df['z'])]
    
    
    def calc_DNA_pos_deviation(self,ref=0,only_average=False):
        """
        Gets the deviation of DNA bp centers from reference
        ref - should be time frame or an a_DNA.df - dataframe
        
        """
#         if hasattr(self, 'contdf_series'):
#             return self.contdf_series
        if isinstance(ref,int):
            ref=self.df_series[self.df_series['Frame']==ref].copy()
            
    
        if(only_average==False):
        #This might probably be speeded up by using scipy.spatial.distance.cdist(dref,d,'sqeuclidean')
        #But currently the timing is not that bad.
            g=self.df_series.groupby(['Frame'])
            dev=g.apply(lambda x: x[['x','y','z']].reset_index()- ref[['x','y','z']].reset_index()).reset_index()
            self.df_series['dx']=dev['x']
            self.df_series['dy']=dev['y']
            self.df_series['dz']=dev['z']
            self.df_series['pos_dev']=np.sqrt(dev['x']**2+dev['y']**2+dev['z']**2)
            print('Added dx,dy,dz and pos_dev columns to self.df_series')


        #Now calculate for average
        self.df['dx']=self.df['x']-ref['x']
        self.df['dy']=self.df['y']-ref['y']
        self.df['dz']=self.df['z']-ref['z']
        self.df['pos_dev']=np.sqrt(self.df['dx']**2+self.df['dy']**2+self.df['dz']**2)

        
        print('Added dx,dy,dz and pos_dev columns to self.df')
        #return self.df_series
        
    
    def calc_DNA_bp_closest(self,ref=0):
        """
        This will for every time frame calculate the closest base pair number in the reference and get its distance.
        """
        d=self.df_series[['x','y','z']].to_numpy()
        d=d.reshape(-1,self.p.DNA_top_length,3)
        if isinstance(ref,int):
            ref=self.df_series[self.df_series['Frame']==ref]
        #else ref should be a df with x,y,z
        dref=ref[['x','y','z']].to_numpy()
        dref=dref.reshape(self.p.DNA_top_length,3)

        mindist=[]
        minbp=[]
        bpshift=[]


        for i in range(0,int(self.df_series.shape[0]/self.p.DNA_top_length)):
            cdist=scipy.spatial.distance.cdist(dref,d[i],'euclidean')
            mindist.extend(np.min(cdist,axis=0))
            m=(cdist == np.min(cdist,axis=0))
            r,c = np.where(m)
            minbp.extend(r+self.p.dyad_top[1]-self.p.DNA_left_length)
            bpshift.extend(r-np.arange(self.p.DNA_top_length))

        self.df_series['BP_mindist']=mindist
        self.df_series['BP_id_closest']=minbp
        self.df_series['BP_id_shift']=bpshift

        print('BP_mindist, BP_id_closest and BP_id_shift added to df_series dataframe')


#         return pm
    
#     def calc_unwrapping(self,threshold=10):
#         """
#         This will 
#         """

    def gen_find_pair_file(self):
        """
        generate find pair file as buffer
        """
        output='input.pdb'+'\n'+'input.out\n'
        output+='    2         # duplex\n%5d         # number of base-pairs\n    1     1    # explicit bp numbering/hetero atoms\n'%len(self.bps)


    #  729  1018   0 #    1 | ....>I:.-72_:[.DA]A-----T[.DT]:..72_:J<....   1.82   1.59  28.50   9.48   4.43
        for n,top_resid in enumerate(self.bps,1):
            bot_resid= self.bp_dict['top'][top_resid]
            resindex_top=self.sel_u.select_atoms('segid %s and resnum %d'%(self.p.DNA_top_segid,top_resid)).residues.resindices[0]+1
            resname_top=self.sel_u.select_atoms('segid %s and resnum %d'%(self.p.DNA_top_segid,top_resid)).residues.resnames[0]
    #        comp_resid=self.p.dyad_bot[1]+self.p.dyad_top[1]-resid
    #             print(resid,comp_resid,self.p.dyad_bot[1],self.p.dyad_top[1])
            resindex_bot=self.sel_u.select_atoms('segid %s and resnum %d'%(self.p.DNA_bot_segid,bot_resid)).residues.resindices[0]+1
            resname_bot=self.sel_u.select_atoms('segid %s and resnum %d'%(self.p.DNA_bot_segid,bot_resid)).residues.resnames[0]
            output+="%5d %5d   0 #%5d | ....>%s:%s_:[%s]%1s-----%1s[%s]:%s_:%s<....   0.00   0.00  00.00   0.00   0.00\n"%(resindex_top,resindex_bot,n,self.p.DNA_top_segid,
                                                                                                                           ('%4d'%top_resid).replace(' ','.'),
                                                                                                                           ('%3s'%resname_top).replace(' ','.'),
                                                                                                                           resname_top[-1],resname_bot[-1],('%3s'%resname_bot).replace(' ','.'),
                                                                                                                           ('%4d'%bot_resid).replace(' ','.'),self.p.DNA_bot_segid)
        #print('DEBUG')
        #print(output)            
        buffer = StringIO()
        buffer.write(output)
        buffer.seek(0)
        return(buffer)
# Old way
#     def gen_find_pair_file(self):
#         """
#         generate find pair file as buffer
#         """
#         output='input.pdb'+'\n'+'input.out\n'
#         output+='    2         # duplex\n%5d         # number of base-pairs\n    1     1    # explicit bp numbering/hetero atoms\n'%len(self.bps)
        
        
# #  729  1018   0 #    1 | ....>I:.-72_:[.DA]A-----T[.DT]:..72_:J<....   1.82   1.59  28.50   9.48   4.43
#         for resid,n in zip(self.bps,range(1,len(self.bps)+1)):
#             resindex_top=self.sel_u.select_atoms('segid %s and resnum %d'%(self.p.DNA_top_segid,resid)).residues.resindices[0]+1
#             resname_top=self.sel_u.select_atoms('segid %s and resnum %d'%(self.p.DNA_top_segid,resid)).residues.resnames[0]
#             comp_resid=self.p.dyad_bot[1]+self.p.dyad_top[1]-resid
# #             print(resid,comp_resid,self.p.dyad_bot[1],self.p.dyad_top[1])
#             resindex_bot=self.sel_u.select_atoms('segid %s and resnum %d'%(self.p.DNA_bot_segid,comp_resid)).residues.resindices[0]+1
#             resname_bot=self.sel_u.select_atoms('segid %s and resnum %d'%(self.p.DNA_bot_segid,comp_resid)).residues.resnames[0]
#             output+="%5d %5d   0 #%5d | ....>%s:%s_:[%s]%1s-----%1s[%s]:%s_:%s<....   0.00   0.00  00.00   0.00   0.00\n"%(resindex_top,resindex_bot,n,self.p.DNA_top_segid,('%4d'%resid).replace(' ','.'),('%3s'%resname_top).replace(' ','.'),resname_top[-1],resname_bot[-1],('%3s'%resname_bot).replace(' ','.'),('%4d'%comp_resid).replace(' ','.'),self.p.DNA_bot_segid)
        
#         print(output)
#         buffer = StringIO()
#         buffer.write(output)
#         buffer.seek(0)
#         return(buffer)

# has to be external to class to be picklable
def X3DNA_analyze_runner(all_vars):
    pdb,frame,sel,bp_ref_file,bps=all_vars
    sel.write(pdb,frames=[frame])
    df=DNAtools.X3DNA_analyze(pdb,ref_fp_id=bp_ref_file)
    df=df.iloc[0:len(bps)]
    df['Frame']=frame
    #df['BPnum_dyad']=bps-p.dyad_top[1]
    #df['segid']=nuclstr_instance.dyad_top[0]
    df['resid']=bps
    return(df) 
    

    
    