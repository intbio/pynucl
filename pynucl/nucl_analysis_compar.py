# -*- coding: utf-8 -*-
"""

Сравнение результатов анализа для разных объектов

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
