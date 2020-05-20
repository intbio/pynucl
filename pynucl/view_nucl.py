import os
import urllib.request
import json
import nglview as nv
import MDAnalysis as mda



def view_nucl(*args,gui=False,chconv={},selection='all',color='bright'):
    '''
    Function shows preview of the nucleosome (chain names alike 1kx5) via NGLview and MDAnalysis
    args   - MDA universe or atom selection, or anything that can be parsed via mda.Universe()
    gui    - shows standard nglview gui (very limited)
    chconv - Dictionary with chain name conformity to map to 1kx5
    selection - mda selection string
    color  - color preset (bright or dull, or dictionary like color={'A':0x94b4d1,'B':0x94d19c,'C':0xd6d989,'D':0xd98989,'DNA':0xd6d6d6})
    
    '''
    ch_conv={'A':'A','B':'B','C':'C','D':'D','E':'E','F':'F','G':'G','H':'H','I':'I','J':'J'}
    ch_conv.update(chconv)
    if(isinstance(args[0],mda.Universe)):
        nuclMD=args[0]
    else:
        nuclMD=mda.Universe(*args)
    #prot = nuclMD.select_atoms("protein")
    sel=nuclMD.select_atoms(selection)
    show=nv.show_mdanalysis(sel,gui=gui)
    if(color=='dull'):
        color={'A':0x94b4d1,'B':0x94d19c,'C':0xd6d989,'D':0xd98989,'DNA':0xd6d6d6}
    elif isinstance(color,dict):
        color=color
    else:
        color={'A':0x020AED,'B':"green",'C':0xE0F705,'D':0xCE0000,'DNA':"grey"}
    
    show.representations = [
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['A']), "color": color['A'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['B']), "color": color['B'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['C']), "color": color['C'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['D']), "color": color['D'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['E']), "color": color['A'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s "%(ch_conv['F']), "color": color['B'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['G']), "color": color['C'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s "%(ch_conv['H']), "color": color['D'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": f"protein and not ( :{ch_conv['A']} :{ch_conv['B']}  :{ch_conv['C']} :{ch_conv['D']} :{ch_conv['E']} :{ch_conv['F']} :{ch_conv['G']} :{ch_conv['H']} )",
        "color":0xEEEEEE ,"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }}, 
    {"type": "cartoon", "params": {
        "sele": "nucleic", "color": color['DNA'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "base", "params": {
        "sele": "nucleic", "color": color['DNA'],
    }},
    ]
    show.camera = 'orthographic'
    return show

def nucl_repr(view,chconv={},selection='all',component=0,color=False):
    ch_conv={'A':'A','B':'B','C':'C','D':'D','E':'E','F':'F','G':'G','H':'H','I':'I','J':'J'}
    ch_conv.update(chconv)
    if(color==False):
        color={'A':0x020AED,'B':"green",'C':0xE0F705,'D':0xCE0000,'DNA':"grey"}
    view.set_representations([
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['A']), "color": color['A'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['B']), "color": color['B'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['C']), "color": color['C'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['D']), "color": color['D'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['E']), "color": color['A'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['F']), "color": color['B'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['G']), "color": color['C'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['H']), "color": color['D'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": "nucleic", "color": color['DNA'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "base", "params": {
        "sele": "nucleic", "color": color['DNA'],
    }},
    ],component=component)
    view.camera = 'orthographic'


