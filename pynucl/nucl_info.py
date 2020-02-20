# -*- coding: utf-8 -*-
"""
parsing nucleosome pdb
анализ файла, определение цепей и их принадлежности к гистоновым вариантам
определение ДНК
визуализация при помощи pytexshade
определение неличия води и ионов
разрешение, b,r фактор, картинка
классификация на варианты и комплексы и структуры из нескольких нуклеосом

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
