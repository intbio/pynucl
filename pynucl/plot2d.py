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

#import matplotlib.pyplot as plt

def plot_coord(inpdata,plane='xy',column='coord'):
    """
    Make plots of a coord data
    """
    
    d=inpdata
    d['x']=d[column].apply(lambda x: x[0])
    d['y']=d[column].apply(lambda x: x[1])
    d['z']=d[column].apply(lambda x: x[2])


    plot=(ggplot(data=d,mapping=aes(x=plane[0], y=plane[1]))
        + geom_point(size=0.1)+xlab('Coordinate X')+ylab('Coordinate Y'))
#     for i in d.segid.unique():
#         plot=plot+geom_path(data=d[])
    g=d.groupby(['segid','Time'])
    for i in g.groups:
        plot=plot+geom_path(data=g.get_group(i))

#         + scale_x_continuous(limits=(0.5,sl+0.5+rof),expand=(0,0.2),name='',breaks=[])
       # + scale_y_continuous(breaks=[0,0.5,1.0])
    
    return plot


# import matplotlib.pyplot as plt
# plt.figure(figsize=(5,5))

# plt.plot(pos[:,0],pos[:,1],label='Dynamics',color='deepskyblue')
# plt.show()


