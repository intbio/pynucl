"""
Contact analysis tools and classes
"""

from pymolint.mol_int import struct2cont # this is for contact analysis
from .nucl_meta_select import create_elem_dict,nucl_sel_expand
#### Analysis classes
class a_contacts(struct2cont):
    """
    Class to analyze contacts in nucleosome
    Get two selections (using expanded meta synataxis)
    """
    def __init__(self,nuclstr_instance,selA,selB=None, **kwargs):
        if 'time' in kwargs:
            tm=kwargs['time']
            newkwargs=kwargs
            del newkwargs['time']
        else:
            tm=nuclstr_instance.time
            newkwargs=kwargs
        super().__init__(nuclstr_instance.u,nucl_sel_expand(selA,nuclstr_instance.nucl_elements),selB=nucl_sel_expand(selB,nuclstr_instance.nucl_elements), time=tm, **newkwargs)
        
        