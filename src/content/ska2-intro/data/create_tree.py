print("> Script for creating & printing a phylogenetic tree!")
print("\t- Checking input file...")
import sys, os
if len(sys.argv) < 2:
    raise RuntimeError("No file provided!")
inf = sys.argv[1]

if not os.path.exists(inf):
    raise RuntimeError("File provided does not exist!")

print("\t- Importing libraries...")
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

print("\t- Reading input...")
aln =  AlignIO.read(inf, "fasta")

print("\t- Calculating distance matrix...")
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

print("\t- Creating tree...")
constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(aln)

print("\t- Drawing tree...")
Phylo.draw_ascii(tree)

print("\nFinished!")
