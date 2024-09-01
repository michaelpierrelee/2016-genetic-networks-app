# coding: utf-8

from . import *

def add_PPI(genome, protein1, protein2):
    """
    add a protein1->protein2 interaction, function adds an edge between nodes "protein1" and "protein2"
    """
    genome.add_edge(protein1, protein2)
    nx.set_edge_attributes(genome, "effect", {(protein1, protein2): "complex"})

def add_PGI(genome, protein, gene, effect, K, hill):
    """
    add a protein->gene interaction (the protein is so a transcription factor)
    """
    genome.add_edge(protein, gene)
    nx.set_edge_attributes(genome, "effect", {(protein, gene): effect})
    nx.set_edge_attributes(genome, "K", {(protein, gene): K})
    nx.set_edge_attributes(genome, "hill", {(protein, gene): hill})

def add_gene(genome, gene, pmax, pmin, degrad, age):
    """
    add a gene producing a protein, don't forget to give a target to this protein
    gene: name of the new gene
    pmax and pmin: production thresholds
    nb = cts["max_nb_of_genes"]
    """
    #add nodes
    gene_id = str(gene).split(".")[1]
    x = 2 + float(gene_id)/10**len(gene_id)
    protein = fckit(x) #index for the node
    genome.add_node(gene)
    genome.add_node(protein)
    #add node attributes
    nx.set_node_attributes(genome, "type", {gene: "gene", protein: "protein"})
    nx.set_node_attributes(genome, "age", {gene: age, protein: age})
    nx.set_node_attributes(genome, "qty", {protein: 0})
    nx.set_node_attributes(genome, "degrad", {protein: degrad})
    nx.set_node_attributes(genome, "pmax", {gene: pmax})
    nx.set_node_attributes(genome, "pmin", {gene: pmin})
    #add edges
    genome.add_edge(gene, protein)
    nx.set_edge_attributes(genome, "effect", {(gene, protein): "transcript"})

def label_complex(protein1, protein2, nb0):
    try: a = str(protein1).split(".")[1]
    except:
        print "ERROR >>>>>>>>>>", protein1
        raise
    b = str(protein2).split(".")[1]
    x = float("3." + nb0 + a + b)
    complex_label = None
    i = 1
    while complex_label != x:
        nb = 10**i
        complex_label = formating(x, nb)
        i += 1
    return complex_label
    
def add_complex(genome, protein1, protein2, coeff1, coeff2, degrad, kon, koff, age):
    """
    add a complex, don't forget to give a target to this complex
    complex_label: name of the new complex
    protein1 and protein2: proteins forming the complex
    a and b: stoechiometric coefficients of the reaction
        a.protein1 + b.protein2 = complex
    """
    #label
    nb0 = ""
    complex_label = label_complex(protein1, protein2, nb0)
    #search if this node already exists
    while genome.has_node(complex_label):
        nb0 += "0"
        complex_label = label_complex(protein1, protein2, nb0)
    #add node
    genome.add_node(complex_label)
    nx.set_node_attributes(genome, "type", {complex_label: "complex"})
    nx.set_node_attributes(genome, "qty", {complex_label: 0})
    nx.set_node_attributes(genome, "degrad", {complex_label: degrad})
    nx.set_node_attributes(genome, "kon", {complex_label: kon})
    nx.set_node_attributes(genome, "koff", {complex_label: koff})
    nx.set_node_attributes(genome, "age", {complex_label: age})
    nx.set_node_attributes(genome, "stoech", { complex_label: {protein1: coeff1, protein2: coeff2} })
    #add edges
    add_PPI(genome, protein1, complex_label)
    add_PPI(genome, protein2, complex_label)
    return complex_label