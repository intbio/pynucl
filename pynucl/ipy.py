from pytexshade import shade
from IPython.display import Image
import tempfile

# class shadedmsa(object):
#     """Class for storing a shaded image"""

#     def __init__(self, msa,shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150,debug=False,startnumber=1):
#         temp = tempfile.NamedTemporaryFile(suffix='.png')
#         if(debug):
#             print("tempfile created: ",temp.name)
#         shade.shade_aln2png(msa, filename=temp.name,shading_modes=shading_modes, features=features,title=title,legend=legend, logo=logo,hideseqs=hideseqs,splitN=splitN,setends=setends,ruler=ruler,show_seq_names=show_seq_names,show_seq_length=show_seq_length,funcgroups=funcgroups,rotate=rotate,threshold=threshold,resperline=resperline,margins=margins, density=density,debug=debug,startnumber=startnumber)
#         self.img=open(temp.name, 'rb').read()

#     def _repr_png_(self):
#         return self.img

