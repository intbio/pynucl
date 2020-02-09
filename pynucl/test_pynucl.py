from pytexshade.shade import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import tempfile
import math
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from io import StringIO
import pkg_resources
import pickle



# def test_shade_aln2png_two_seqs():
# 	os.system('mkdir -p test_results')
# 	print("Testing a small two sequence alignment")
# 	human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
# 	xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
# 	# human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
# 	msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2AZ',name='H2AZ')])
# 	# features=get_hist_ss_in_aln_for_shade(msa,below=True)
# 	features=[{'style':'fill:$\\uparrow$','sel':[5,10],'text':'test'}]
# 	print(features)
# 	shade_aln2png(msa,filename='test_results/2seqs.png',shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
# 	size=os.path.getsize('test_results/2seqs.png')
# 	assert size>10000, "output png filesize too small, looks that nothing was produced"

# def test_5z3l():
# 	DATA_PATH = pkg_resources.resource_filename('pytexshade', 'data/')

# # 	for dirname, dirnames, filenames in os.walk('.'):
# # # print path to all subdirectories first.
# # 		for subdirname in dirnames:
# # 			print(os.path.join(dirname, subdirname))
# # # print path to all filenames.
# # 		for filename in filenames:
# # 			print(os.path.join(dirname, filename))
# # 	print(os.getcwd())
# # 	print(os.path.dirname(os.path.realpath(__file__)))
# 	os.system('mkdir -p test_results')
# 	msa_dict = pickle.load( open(os.path.join(DATA_PATH,"5z3l_msa_dict.p"), "rb" ) )
# 	features_dict = pickle.load( open(os.path.join(DATA_PATH,"5z3l_features_dict.p"), "rb" ) )

# 	for s in 'ABCDEFGHO':
# 		msa=msa_dict[s]
# 		features=features_dict[s]
# 		print("Testing chain %s in 5z3l test"%s)
# 		shade_aln2png(msa,\
#  filename='test_results/5z3l_%s.png'%s,\
#  shading_modes=['similar'],# list of shading modes according to TexShade, supported are "similar", "hydropathy_functional", "chemical_functional", "structure_functional", "charge_functional", "diverse"\
#  legend=False,features=features,title='',\
#  logo=False, #SeqLogo \
#  hideseqs=False,#activate \hidseqs command \ 
#  splitN=20, #alignment will be split into splitN batches\
#  setends=[], # a section of alignment will be displayed between setends[0] and setends[1]\
#  ruler=True, # Add a ruler\
#  show_seq_names=True,# Show sequence names\
#  funcgroups=None, # Hilight functional groups, Tex code should be inserted example funcgroups="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{yyy}{GA}{White}{Blue}{upper}{up}" \
#  show_seq_length=False #Show sequence length \
# )
# 		size=os.path.getsize('test_results/5z3l_%s.png'%s)
# 		assert size>10000, "output png filesize too small, looks that nothing was produced"


# def test_logo():
# 	#Prepare a Multiple Sequence Alignment in biopython
# 	human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
# 	xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
# 	msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2AZ',name='H2AZ')])
# #We need to describe feautures
# #See also TexShade docs http://mirrors.mi.ras.ru/CTAN/macros/latex/contrib/texshade/texshade.pdf
# 	features=[
#     {'style':'fill:$\\uparrow$', #Styles - all TexShade feature types 'helix','loop','-->','---','<--',',-,' AND also "frameblock", "shaderegion" or "shadeblock", see TexShade docs for more info\
#      #'position':'top',# Position of feature annotation also 'bottom','ttop','bbottom', etc. if no - automatic'\
#      #'seqref':1, #number of sequence for selections - default consensus\
#      'sel':[5,10], # selection with 0-based numbering, both start and end are included!\
#      #'color':'Red',# color for for frame block\
#      #'thickness':1.5 #-1.5pt also for frameblock\
#      #'rescol':'Black',  #this is residue color for shaderegion style, see TexShade docs \shaderegion{ seqref  }{ selection }{ res.col. }{ shad.col. }\
#      #'shadcol':'Green',#this is shding color for shaderegion style, see TexShade docs \shaderegion{ seqref  }{ selection }{ res.col. }{ shad.col. }\
#      #funcgroup example fg="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{xxx}{GA}{White}{Blue}{upper}{up}"\
#      'text':'test'},
#      {'style':'helix', \
#      'sel':[25,30], # selection with 0-based numbering, both start and end are included!\
#      'text':'$\\alpha 1$'},
#      {'style':'-->', \
#      'sel':[35,40], # selection with 0-based numbering, both start and end are included!\
#      'text':'$\\beta 2$'},
# ]
# #Shade
# 	shade_aln2png(msa,\
#  	filename='test_results/logo.png',\
#  	 debug=False,\
#      shading_modes=['similar'],# list of shading modes according to TexShade, supported are "similar", "hydropathy_functional", "chemical_functional", "structure_functional", "charge_functional", "diverse"\
#      legend=False,features=features,title='',\
#      logo=True, #SeqLogo \
#      hideseqs=False,#activate \hidseqs command \ 
#      splitN=20, #alignment will be split into splitN batches\
#      setends=[], # a section of alignment will be displayed between setends[0] and setends[1]\
#      ruler=False, # Add a ruler\
#      show_seq_names=True,# Show sequence names\
#      funcgroups=None, # Hilight functional groups, Tex code should be inserted example funcgroups="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{yyy}{GA}{White}{Blue}{upper}{up}" \
#      show_seq_length=True #Show sequence length \
#     )

# 	size=os.path.getsize('test_results/logo.png')
# 	assert size>10000, "output png filesize too small, looks that nothing was produced"


# def test_long_msa():
# 	from Bio import AlignIO
# 	from urllib.request import urlopen 
# 	import tempfile
# 	temp = tempfile.NamedTemporaryFile()
# 	pfam_path="https://pfam.xfam.org/family/PF00125/alignment/seed"
# 	data = urlopen(pfam_path).read()
# 	# print(data)
# 	with open(temp.name,'wb') as f:
# 		f.write(data)
# 	alignment = AlignIO.read(temp.name, "stockholm")

# #Shade
# 	shade_aln2png(alignment,\
# 	filename='test_results/large.png',\
#      shading_modes=['similar'],# list of shading modes according to TexShade, supported are "similar", "hydropathy_functional", "chemical_functional", "structure_functional", "charge_functional", "diverse"\
#      legend=False,features=[],title='',\
#      logo=True, #SeqLogo \
#      hideseqs=False,#activate \hidseqs command \ 
#      splitN=20, #alignment will be split into splitN batches\
#      setends=[], # a section of alignment will be displayed between setends[0] and setends[1]\
#      ruler=False, # Add a ruler\
#      show_seq_names=True,# Show sequence names\
#      funcgroups=None, # Hilight functional groups, Tex code should be inserted example funcgroups="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{yyy}{GA}{White}{Blue}{upper}{up}" \
#      show_seq_length=True #Show sequence length \
#     )

# 	size=os.path.getsize('test_results/large.png')
# 	assert size>10000, "output png filesize too small, looks that nothing was produced"


