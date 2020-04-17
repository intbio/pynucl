import os
import urllib.request
import json
import nglview as nv
import MDAnalysis as mda


def view_nucl(*args,gui=False,chconv={},selection='all',color=False):
    ch_conv={'A':'A','B':'B','C':'C','D':'D','E':'E','F':'F','G':'G','H':'H','I':'I','J':'J'}
    ch_conv.update(chconv)
    if(isinstance(args[0],mda.Universe)):
        nuclMD=args[0]
    else:
        nuclMD=mda.Universe(*args)
    #prot = nuclMD.select_atoms("protein")
    sel=nuclMD.select_atoms(selection)
    show=nv.show_mdanalysis(sel,gui=gui)
    if(color==False):
        color={'A':0x020AED,'B':"green",'C':0xE0F705,'D':0xCE0000,'DNA':"grey"}
    show.representations = [
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
        "sele": ":%s "%(ch_conv['F']), "color": color['B'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['G']), "color": color['C'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s "%(ch_conv['H']), "color": color['D'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
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


