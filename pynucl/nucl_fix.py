# -*- coding: utf-8 -*-
"""
rename/renamber histone/dna
orient NRF
output pdb and ful sequence

"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import tempfile
import math
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from io import StringIO
