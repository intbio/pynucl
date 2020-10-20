import os
import urllib.request
import json
import nglview as nv
import MDAnalysis as mda


def view_nucl(*args,legend=True,gui=False,onlyNcp=False,chconv={},selection='all',color='bright'):
    '''
    Function shows preview of the nucleosome (chain names alike 1kx5) via NGLview and MDAnalysis
    args   - MDA universe or atom selection, or anything that can be parsed via mda.Universe()
    legend - show legend
    gui    - shows standard nglview gui (very limited)
    onlyNcp - show only ncp
    chconv - Dictionary with chain name conformity to map to 1kx5
    selection - mda selection string
    color  - color preset (bright or dull, or dictionary like color={'H3':"#94b4d1",'H4':"#94d19c",'H2A':"#d6d989",'H2B':"#d98989",'DNA':"#d6d6d6"})
    
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
        color={'H3':"#94b4d1",'H4':"#94d19c",'H2A':"#d6d989",'H2B':"#d98989",'DNA':"#d6d6d6"}
    elif isinstance(color,dict):
        color=color
    else:
        color={'H3':"#020AED",'H4':"#008000",'H2A':"#E0F705",'H2B':"#CE0000",'DNA':"#808080"}
    
    
    representations=[
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['A']), "color": color['H3'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['B']), "color": color['H4'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['C']), "color": color['H2A'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['D']), "color": color['H2B'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['E']), "color": color['H3'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s "%(ch_conv['F']), "color": color['H4'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s"%(ch_conv['G']), "color": color['H2A'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }},
    {"type": "cartoon", "params": {
        "sele": ":%s "%(ch_conv['H']), "color": color['H2B'],"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
    }}, 
    {"type": "cartoon", "params": {
        "sele": "nucleic", "color": color['DNA'],"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "base", "params": {
        "sele": "nucleic", "color": color['DNA'],
    }},
    ]
    if not onlyNcp:
        representations.append({"type": "cartoon", "params": {
                                "sele": f"protein and not ( :{ch_conv['A']} :{ch_conv['B']}  :{ch_conv['C']} :{ch_conv['D']} :{ch_conv['E']} :{ch_conv['F']} :{ch_conv['G']} :{ch_conv['H']} )",
                                "color":0xEEEEEE ,"aspectRatio":2, 'radiusScale':4,'radiusType':'sstruc',"capped":True,'subdiv':10,'diffuseInterior':False,'useInteriorColor':False
                                }})
    show.representations = representations
    show.camera = 'orthographic'
    def add_mda_representaton(representation,sel,name=None,color='red',**kwargs):
        """
        Adds ngl representation 
        see http://nglviewer.org/ngl/api/manual/molecular-representations.html
        from mda selection string
        all ngl kwargs can be passed!
        """
        #indeces are numbered the same way as in nglview
        atom_ids=nuclMD.select_atoms(sel).indices
        show.add_representation(representation,sele='@'+','.join(atom_ids.astype(str)),color=color,**kwargs)
        if isinstance(name,str):
            anotate_line(name,show.current_line,color)
            show.current_line+=1
    show.add_mda_representaton=add_mda_representaton
    
    def anotate_line(text,line,color,fontsize=20):
        code=f"""
        var $text = $("<div></div>")
                    .css("position","absolute")
                    .css("top","{fontsize*line+1}px")
                    .css("left","1px")
                    .css("padding","2px 5px 2px 5px")
                    .css("opacity","1.0")
                    .css("font-size","{fontsize}px")
                    .css("color","{color}")
                    .css("text-shadow","-1px 0 black, 0 1px black, 1px 0 black, 0 -1px black")
                    .css("font-family","verdana")
                    .css("font-weight","bold")
                    .appendTo(this.$container);
        $text.text("{text}")
        """
        show._execute_js_code(code)
    show.current_line=0
    if legend:
        for text,color in color.items():
            anotate_line(text,show.current_line,color)
            show.current_line+=1    
    def add_legend_line(text,color):
        """
        adds new line to legend
        text
        color
        """
        anotate_line(text,show.current_line,color)
        show.current_line+=1
    show.add_legend_line=add_legend_line
    
    def set_size(x,y):
        show._remote_call('setSize',target='Widget',args=[f'{x}px',f'{y}px'])
        show.center()
    show.set_size=set_size
    
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


