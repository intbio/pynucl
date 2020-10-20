# Some visualization helper classes or functions here
# including wrappers for seqplot and pytexshade

from plotnine.utils import resolution
from plotnine.doctools import document
from plotnine.geoms.geom_tile import geom_tile
from plotnine import aes, theme

from pytexshade import ipyshade,shade
from plotnine import ggplot,geom_rect, geom_point, aes, stat_smooth,geom_bar, xlim, ylim, facet_wrap, theme_bw,theme_xkcd, geom_line, geom_tile
from plotnine import facet_wrap, theme, scale_y_continuous,scale_x_continuous, theme_bw,theme_classic, theme_dark, theme_light, theme_matplotlib, theme_minimal, theme_seaborn, theme_void, geom_rect, element_text
from plotnine.labels import ggtitle
from plotnine.scales import scale_fill_manual
from plotnine.guides import guides
from plotnine.guides import guide_legend

from seqplot.p9tools import geom_seq_x
from seqplot.mpltools import add_seq_x



import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio import SeqIO
from Bio import Entrez
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna, generic_protein

import logging
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)s %(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

def multi_segid(func):
    def for_every_segid(data,*args, **kwargs):
        segids=data.segid.unique()
        if len(segids)==1:
            return func(data,*args, **kwargs)
        if len(segids)>1:
            retobj=[]
            for s in segids:
                retobj.append(func(data[data.segid.eq(s)],*args,**kwargs))
        return retobj
    return for_every_segid

@multi_segid
def plot_bar(data,nuclstr,column='value',factor=None,ymin=None,ymax=None,stat='identity',dpi=300,features=None,feature_types=['all'],add_features=[],funcgroups=None,shading_modes=['charge_functional'],usd=False,right_overhang_fix=None,debug=False,startnumber=1,cropseq=(0,None),aspect_ratio=None,reverse_seq=False,double_seq=False,transparent=True,fill_params=None,bar_position='stack',title=None):
    """
    A wrapper function to make a plot of data with bars along the sequnce
    input should be a dataframe with resid, segid column and 'value' 
    This one is inspired by seqplot/seqplot/pdb_plot.py
    """
    
    segid=data['segid'].values[0]
    
    if title is None:
        title="Segid: %s, Type: %s"%(segid,nuclstr.components[segid]['type'])
    
    seq=Seq(str(nuclstr.seqs[segid]['fullseq']),generic_protein \
                if nuclstr.components[segid]['entity'] is 'DNA' or 'histone' or 'protein' else generic_dna)
    msar=MultipleSeqAlignment([SeqRecord(seq=seq,id=nuclstr.components[segid]['type']+':'+segid,\
                                         name=nuclstr.components[segid]['type']+':'+segid)])
    if(reverse_seq):
        logger.info("Experimental feature will reverse the sequence")
        msar[0].seq=msar[0].seq[::-1]

    if double_seq:
          msar.add_sequence('reverse',str(msar[0].seq[::-1]))

        
    msar=msar[:,cropseq[0]:cropseq[1]]
        
    
#     print("Seq to plot:",msar)
             
    #We need to get starting residue, currently for DNA chains only cifseq gets it correctly
    resid_start=nuclstr.seqs[segid]['resid_start']
    
    logger.debug("Starting resid",resid_start)
    

    overhang=nuclstr.seqs[segid]['overhangL']
    
    datafixed=data.copy()
    datafixed.loc[:,'resid']=datafixed.loc[:,'resid']-resid_start+overhang+1-cropseq[0]

    
    sl=len(msar[0].seq)

#     fn=shade.seqfeat2shadefeat(msar,feature_types=feature_types,force_feature_pos='bottom',debug=debug)
    if features is None:
        fn=nuclstr.shading_features[segid]
    else:
        fn=features
    fn2=[]
    for i in fn:
        if (i['style'] in feature_types) or ('all' in feature_types) :
            fn2.append(i)
            
    fn2.extend(add_features)
    if usd:
        ruler='top'
    else:
        ruler=None
    shaded=ipyshade.shadedmsa4plot(msar,features=fn2,shading_modes=shading_modes,debug=debug,startnumber=startnumber,setends=[startnumber-2,sl+startnumber+2],funcgroups=funcgroups,ruler=ruler,density=200)
        
    #If sl%10=10 se will have a ruler number hanging beyond the sequence image, and we need to correct for that.
    if right_overhang_fix is None:
        if sl%10==0:
            if sl<100:
                rof= 0.1
            else:
                rof=0.5
        else:
            rof=0
    else:
        rof=right_overhang_fix
    if (not aspect_ratio is None ):
        ar=aspect_ratio
    else:
        ar=0.2*100./sl
#     print(datafixed)
    plot=(ggplot(data=datafixed,mapping=aes(x='resid', y=column))
#         + geom_point(size=0.1)
#           +geom_bar(stat='identity',width=0.5,mapping=aes(fill=factor))
        + scale_x_continuous(limits=(0.5,sl+0.5+rof),expand=(0,0.2),name='',breaks=[])
       # + scale_y_continuous(breaks=[0,0.5,1.0])
        + theme_light()+theme(aspect_ratio=ar,dpi=dpi,plot_margin=0,text=element_text(size=6), legend_key_size=5 ,legend_position='bottom',legend_direction='horizontal'))
    #+ facet_wrap('~ segid',dir='v') +guides(color=guide_legend(ncol=10))
    if factor is None:
        plot=plot+geom_bar(stat=stat,width=0.5)
    else:
        plot=plot+geom_bar(stat=stat,width=0.5,mapping=aes(fill=factor),position=bar_position)
        
    if fill_params is not None:
        plot=plot+scale_fill_manual(**fill_params)
    
    if not usd:
        if (ymax is not None) :
            plot=plot+scale_y_continuous(limits=(None,ymax))
    else:
        if (ymin is not None) :
            plot=plot+scale_y_continuous(limits=(ymin,None))
    
    if ymax is None:
        ymax=data[column].max()
    if ymin is None:
        ymin=data[column].min()
#     print(ymax)
    plot = plot + geom_seq_x(seqimg=shaded.img,\
                   xlim=(1,sl+rof),ylim=(ymin,ymax),usd=usd,aspect_ratio=ar,transparent=transparent)+ggtitle(title)
    
    
    return plot
        
    

@multi_segid    
def plot_line(data,nuclstr,columns=['value'],ymin=None,ymax=None,dpi=300,features=None,feature_types=['all'],add_features=[],funcgroups=None,shading_modes=['charge_functional'],right_overhang_fix=None,debug=False,startnumber=1,cropseq=(0,None),aspect_ratio=None,reverse_seq=False,transparent=True,xshift=0):
    """
    A wrapper function to make a plot of data with bars along the sequnce
    input should be a dataframe with resid, segid column and 'value' 
    This one is inspired by seqplot/seqplot/pdb_plot.py
    funcgroup example fg="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{xxx}{GA}{White}{Blue}{upper}{up}"
    """
    if isinstance(columns,str):
        columns=[columns]
    segid=data['segid'].values[0]
    
    title="Segid: %s, Type: %s"%(segid,nuclstr.components[segid]['type'])

    seq=Seq(str(nuclstr.seqs[segid]['fullseq']),generic_protein \
                if nuclstr.components[segid]['entity'] is 'DNA' or 'histone' or 'protein' else generic_dna)
    msar=MultipleSeqAlignment([SeqRecord(seq=seq,id=nuclstr.components[segid]['type']+':'+segid,\
                                         name=nuclstr.components[segid]['type']+':'+segid)])
    if(reverse_seq):
        logger.info("Experimental feature will reverse the sequence")
        msar[0].seq=msar[0].seq[::-1]
        
    msar=msar[:,cropseq[0]:cropseq[1]]

    
#     print("Seq to plot:",msar)
             
    #We need to get starting residue, currently for DNA chains only cifseq gets it correctly
    resid_start=nuclstr.seqs[segid]['resid_start']
    
    logger.debug("Starting resid %d"%int(resid_start))
    

    overhang=nuclstr.seqs[segid]['overhangL']
    
    datafixed=data.copy()
    datafixed.loc[:,'resid']=datafixed.loc[:,'resid']-resid_start+overhang+1-cropseq[0]+xshift

#     print(datafixed)
    sl=len(msar[0].seq)

#     fn=shade.seqfeat2shadefeat(msar,feature_types=feature_types,force_feature_pos='bottom',debug=debug)
    if features is None:
        fn=nuclstr.shading_features[segid]
    else:
        fn=features
    fn2=[]
    for i in fn:
        if (i['style'] in feature_types) or ('all' in feature_types) :
            fn2.append(i)
            
    fn2.extend(add_features)
    shaded=ipyshade.shadedmsa4plot(msar,features=fn2,shading_modes=shading_modes,debug=debug,startnumber=startnumber,setends=[startnumber-2,sl+startnumber+2],funcgroups=funcgroups,density=200)
        
    #If sl%10=10 se will have a ruler number hanging beyond the sequence image, and we need to correct for that.
    if right_overhang_fix is None:
        if sl%10==0:
            if sl<100:
                rof= 0.1
            else:
                rof=0.5
        else:
            rof=0
    else:
        rof=right_overhang_fix
    if (not aspect_ratio is None ):
        ar=aspect_ratio
    else:
        ar=0.15*100./sl
        
    md=pd.melt(datafixed,id_vars=['segid','resid'],value_vars=columns)
#     print(md)
#     print(md)
#     print(md['variable'])
    plot=(ggplot(data=md,mapping=aes(x='resid', y='value'))
        + geom_point(aes(color='variable'),size=0.1)+geom_line(aes(color='variable'),stat='identity')
        + scale_x_continuous(limits=(0.5,sl+0.5+rof),expand=(0,0.2),name='',breaks=[])
#         + scale_y_continuous()
        + theme_light()+theme(aspect_ratio=ar,dpi=dpi,plot_margin=0)) #+ facet_wrap('~ segid',dir='v')

    if ymax is not None:
        plot=plot+scale_y_continuous(limits=(None,ymax))
    
    if ymin is None:
        ymin=md['value'].min()
    if ymax is None:
        ymax=md['value'].max()
    plot = plot + geom_seq_x(seqimg=shaded.img,\
                   xlim=(1,sl+rof),ylim=(ymin,ymax),aspect_ratio=ar,transparent=transparent)+ggtitle(title)
    

    
    return plot



@multi_segid    
def plot_line_mpl(data,nuclstr,column='value',error=None,ref=0,ymin=None,ymax=None,dpi=300,features=None,feature_types=['all'],add_features=[],funcgroups=None,shading_modes=['charge_functional'],right_overhang_fix=None,debug=False,startnumber=1,cropseq=(0,None),aspect_ratio=None,reverse_seq=False,transparent=True,xshift=0,title=None,color='blue',color_ref='red',figsize=(15,3),ax=None):
    """
    A wrapper function to make a plot of data with bars along the sequnce
    input should be a dataframe with resid, segid column and 'value' 
    This one is inspired by seqplot/seqplot/pdb_plot.py
    funcgroup example fg="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{xxx}{GA}{White}{Blue}{upper}{up}"
    """
    if ax is None:
        fig,ax=plt.subplots(figsize=figsize,dpi=dpi) 
    segid=data['segid'].values[0]
    
   
    ax.get_figure().suptitle(title)
    
    title="Segid: %s, Type: %s"%(segid,nuclstr.components[segid]['type'])

    seq=Seq(str(nuclstr.seqs[segid]['fullseq']),generic_protein \
                if nuclstr.components[segid]['entity'] is 'DNA' or 'histone' or 'protein' else generic_dna)
    msar=MultipleSeqAlignment([SeqRecord(seq=seq,id=nuclstr.components[segid]['type']+':'+segid,\
                                         name=nuclstr.components[segid]['type']+':'+segid)])
    if(reverse_seq):
        logger.info("Experimental feature will reverse the sequence")
        msar[0].seq=msar[0].seq[::-1]
        
    msar=msar[:,cropseq[0]:cropseq[1]]

    
#     print("Seq to plot:",msar)
             
    #We need to get starting residue, currently for DNA chains only cifseq gets it correctly
    resid_start=nuclstr.seqs[segid]['resid_start']
    
    logger.debug("Starting resid %d"%int(resid_start))
    

    overhang=nuclstr.seqs[segid]['overhangL']
    
    datafixed=data.copy()
    datafixed.loc[:,'resid']=datafixed.loc[:,'resid']-resid_start+overhang+1-cropseq[0]+xshift

    if 'Frame' in data.columns:
        g=datafixed.groupby(['segid','Frame'])
    else:
        g=datafixed.groupby(['segid'])   
    
#     print(data.datafixed)


#     print(datafixed)
    sl=len(msar[0].seq)

#     fn=shade.seqfeat2shadefeat(msar,feature_types=feature_types,force_feature_pos='bottom',debug=debug)
    if features is None:
        fn=nuclstr.shading_features[segid]
    else:
        fn=features
    fn2=[]
    for i in fn:
        if (i['style'] in feature_types) or ('all' in feature_types) :
            fn2.append(i)
            
    fn2.extend(add_features)
    shaded=ipyshade.shadedmsa4plot(msar,features=fn2,shading_modes=shading_modes,debug=debug,startnumber=startnumber,setends=[startnumber-2,sl+startnumber+2],funcgroups=funcgroups,density=200)
        
    #If sl%10=10 se will have a ruler number hanging beyond the sequence image, and we need to correct for that.
    if right_overhang_fix is None:
        if sl%10==0:
            if sl<100:
                rof= 0.1
            else:
                rof=0.5
        else:
            rof=0
    else:
        rof=right_overhang_fix
    if (not aspect_ratio is None ):
        ar=aspect_ratio
    else:
        ar=0.15*100./sl
    
        
    for i in g.groups:
        ax.plot(g.get_group(i)['resid'].values,g.get_group(i)[column].values,color=color)

#     md=pd.melt(datafixed,id_vars=['segid','resid'],value_vars=columns)
#     plot=(ggplot(data=md,mapping=aes(x='resid', y='value'))
#         + geom_point(aes(color='variable'),size=0.1)+geom_line(aes(color='variable'),stat='identity')
#         + scale_x_continuous(limits=(0.5,sl+0.5+rof),expand=(0,0.2),name='',breaks=[])
#         + theme_light()+theme(aspect_ratio=ar,dpi=dpi,plot_margin=0)) #+ facet_wrap('~ segid',dir='v')

#     if ymin is None:
#         ymin=md['value'].min()
#     plot = plot + geom_seq_x(seqimg=shaded.img,\
#                    xlim=(1,sl+rof),ylim=(ymin,md['value'].max()),aspect_ratio=ar,transparent=transparent)+ggtitle(title)
    
    if error is not None:
        for i in g.groups:
            ax.errorbar(g.get_group(i)['resid'].values,g.get_group(i)[column].values,g.get_group(i)    [error].values,color=color)
    
    if isinstance(ref,int):
        r=data[data['Frame']==ref]
    else:
        r=ref
    if r is not None:
        print("Plotting ref")
        rfixed=r.copy()
        rfixed.loc[:,'resid']=rfixed.loc[:,'resid']-resid_start+overhang+1-cropseq[0]+xshift
        g=rfixed.groupby(['segid'])
        for i in g.groups:
            ax.plot(g.get_group(i)['resid'].values,g.get_group(i)[column].values,color=color_ref)

    if ymin is None:
        if ref is not None:
            ymin=min(datafixed[column].min(),rfixed[column].min())
        else:
            ymin=datafixed[column].min()
    if ymax is None:
        if ref is not None:
            ymax=max(datafixed[column].max(),rfixed[column].max(),)
        else:
            ymax=datafixed[column].max()

#     print(ymin,datafixed[column].max())
    ax=add_seq_x(ax,shaded.img,xlim=(1,sl+rof),ylim=(ymin,ymax))
    return ax