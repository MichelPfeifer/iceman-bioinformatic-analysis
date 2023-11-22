from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib.pyplot as plt

# Code for phylogenetic tree construction derived from: https://github.com/taylor-lindsay/phylogenetics
align = AlignIO.read("phylo/sequences/alignment.fasta", "fasta")

calculator = DistanceCalculator("identity")
dist = calculator.get_distance(align)

constructor = DistanceTreeConstructor()
tree = constructor.upgma(dist)

# function get_label() derived from: https://stackoverflow.com/questions/33147651/biopython-phylogenetic-tree-edit-labels-in-svg-file
def get_label(leaf):
    if not ((leaf.name).startswith("Inner")):
        return leaf.name

def construct_tree(tree, out):
    fig = plt.figure(figsize=(80, 40), dpi=150)
    # Code derived from: https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
    plt.rcParams.update({"font.size": 40})
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False, label_func=get_label)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return

construct_tree(tree, "phylo/results/phylo_tree.png")
