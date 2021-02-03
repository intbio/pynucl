# Some visualization helper classes or functions here
# for polotting of 2D data

from plotnine.utils import resolution
from plotnine.doctools import document
from plotnine.geoms.geom_tile import geom_tile
from plotnine import aes, theme

from pytexshade import ipyshade,shade
from plotnine import ggplot,geom_rect, geom_point, aes, stat_smooth,geom_bar, xlim, ylim, facet_wrap, theme_bw,theme_xkcd, geom_line, geom_tile
from plotnine import facet_wrap, theme, scale_y_continuous,scale_x_continuous, theme_bw,theme_classic, theme_dark, theme_light, theme_matplotlib, theme_minimal, theme_seaborn, theme_void, geom_rect,xlab,ylab,geom_path
from plotnine.labels import ggtitle

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd


#import matplotlib.pyplot as plt

def plot_coord_gg(inpdata,plane='xy',column='coord'):
    """
    Make plots of a coord data
    """
    
    d=inpdata
    d['x']=d[column].apply(lambda x: x[0])
    d['y']=d[column].apply(lambda x: x[1])
    d['z']=d[column].apply(lambda x: x[2])


    plot=(ggplot(data=d,mapping=aes(x=plane[0], y=plane[1]))
        + geom_point(size=0.1)+xlab('Coordinate X')+ylab('Coordinate Y'))
    g=d.groupby(['segid','Frame'])
    for i in g.groups:
        plot=plot+geom_path(data=g.get_group(i))

    
    return plot


def plot_coord(inpdata,plane='xy',column='coord',ref=0,figsize=(5,5),ax=None,color='blue',color_ref='red'):
    """
    Make plots of a coord data
    """
    c={'x':0,'y':1,'z':2}
    if ax is None:
        fig,ax=plt.subplots(figsize=figsize)
    if 'Frame' in inpdata.columns:
        g=inpdata.groupby(['segid','Frame'])
    else:
        g=inpdata.groupby(['segid'])        
    for i in g.groups:
        ax.plot(g.get_group(i)['coord'].apply(lambda x: x[c[plane[0]]]).values,g.get_group(i)['coord'].apply(lambda x: x[c[plane[1]]]).values,color=color)
    
    if isinstance(ref,int):
        r=inpdata[inpdata['Frame']==ref]
    else:
        r=ref
    if r is not None:
        g=r.groupby(['segid'])
        for i in g.groups:
            ax.plot(g.get_group(i)['coord'].apply(lambda x: x[c[plane[0]]]).values,g.get_group(i)['coord'].apply(lambda x: x[c[plane[1]]]).values,color=color_ref)
     
    v1=inpdata['coord'].apply(lambda x: x[c[plane[0]]]).values
    v2=inpdata['coord'].apply(lambda x: x[c[plane[1]]]).values
    sz=np.max([np.max(v1)-np.min(v1),np.max(v2)-np.min(v2)])
    ax.set_xlim((np.max(v1)+np.min(v1))/2-sz/2*1.05,(np.max(v1)+np.min(v1))/2+sz/2*1.05)
    ax.set_ylim((np.max(v2)+np.min(v2))/2-sz/2*1.05,(np.max(v2)+np.min(v2))/2+sz/2*1.05)
    ax.set_xlabel('Coordinate %s, A'%plane.upper()[0])
    ax.set_ylabel('Coordinate %s, A'%plane.upper()[1])

    return ax
    
#     plt.show()

def plot_coord_fast(inpdata,plane='xy',column='coord',ref=0,figsize=(5,5),ax=None,color='blue',color_ref='red',alpha=1.0):
    """
    Make plots of a coord data, faster way using line collections and vectorized pandas-numpy operations
    TODO: add optional animation output http://louistiao.me/posts/notebooks/save-matplotlib-animations-as-gifs/
    """
    c={'x':0,'y':1,'z':2}
    if ax is None:
        fig,ax=plt.subplots(figsize=figsize)
    if 'Frame' in inpdata.columns:
        g=inpdata.groupby(['segid','Frame'])
        array=np.array(list(g['coord'].apply(pd.DataFrame.to_numpy).apply(np.vstack).to_numpy()))
        lines=array[:,:,[c[plane[0]],c[plane[1]]]]
        ln_coll=LineCollection(lines,color=color,alpha=alpha)
        ax.add_collection(ln_coll)
    else:
        g=inpdata.groupby(['segid'])
        for i in g.groups:
            coords=np.stack(g.get_group(i)['coord'].values,axis=0).T
            ax.plot(coords[c[plane[0]]],coords[c[plane[1]]],color=color,alpha=alpha)

    #for i in g.groups:
    #    coords=np.stack(g.get_group(i)['coord'].values,axis=0).T
    #    ax.plot(coords[c[plane[0]]],coords[c[plane[1]]],color=color)
    
    if isinstance(ref,int):
        r=inpdata[inpdata['Frame']==ref]
    else:
        r=ref
    if r is not None:
        g=r.groupby(['segid'])
        for i in g.groups:
            coords=np.stack(g.get_group(i)['coord'].values,axis=0).T
            ax.plot(coords[c[plane[0]]],coords[c[plane[1]]],color=color_ref)
     
    v1=inpdata['coord'].apply(lambda x: x[c[plane[0]]]).values
    v2=inpdata['coord'].apply(lambda x: x[c[plane[1]]]).values
    sz=np.max([np.max(v1)-np.min(v1),np.max(v2)-np.min(v2)])
    ax.set_xlim((np.max(v1)+np.min(v1))/2-sz/2*1.05,(np.max(v1)+np.min(v1))/2+sz/2*1.05)
    ax.set_ylim((np.max(v2)+np.min(v2))/2-sz/2*1.05,(np.max(v2)+np.min(v2))/2+sz/2*1.05)
    ax.set_xlabel('Coordinate %s, A'%plane.upper()[0])
    ax.set_ylabel('Coordinate %s, A'%plane.upper()[1])

    return ax

