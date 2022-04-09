import sys
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import *
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp
from Bio.SeqUtils.ProtParam import ProteinAnalysis

if len(sys.argv) < 2:
    print("Not enough parameters.")
    exit(64)

pLen = int(sys.argv[1])
cDNA = Seq(sys.argv[2])


def get_data(seq):
    mt = MeltingTemp.Tm_Wallace(seq)
    gc = GC(seq)
    mw = molecular_weight(seq)
    ssf = ProteinAnalysis(str(seq)).secondary_structure_fraction()
    return [str(seq), mt, gc, mw, ssf[0], ssf[1], ssf[2]]


# Primer 1
cmp = cDNA.complement()
p1 = cmp[0:pLen]
p1data = get_data(p1)

# Primer 2
rcs = cDNA.reverse_complement()
p2 = rcs[0:pLen]
p2data = get_data(p2)

# Compile into DataFrame
data = {'Primer 1': p1data, 'Primer 2': p2data}
df = pd.DataFrame(data, index=['Seq', 'Tm', 'Mol. Weight', 'GC%', 'fHelix', 'fTurn', 'fSheet'])

# Print
print(df)
