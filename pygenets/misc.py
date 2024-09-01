# coding: utf-8

from . import *

def set_inputs(genome, inputs):
    """
	change the values of inputs of genome
    inputs = {"Input1": value, "Input2": value...}
    """
    for k in inputs.keys():
        nx.set_node_attributes(genome, "input_qty", {
                k: inputs[k]
            })

def update_age(genome):
    """add 1 to each node of genome"""
    for n, dic in genome.nodes(data = True):
        nx.set_node_attributes(genome, "age", {n: dic["age"] + 1})
            
def toss():
    if random.randint(0, 1) == 0: return True
    else: return False
    
def proba_uni(a,b,d):
    """
    uniform probability with a and b the extremities of the invertal
    and d the length of the interval
        a, b and d = integers
    return True if the toss is in the subset, False if not
    """
    if random.randint(a, b) <= d: return True
    else: return False
    
def gene_name(genome):
    #take id of genes and sort the list
    nodes = [int(str(node).split(".")[1]) for (node, dic) in genome.nodes(data = True) if dic["type"] == "gene"]
    if nodes == []: nodes = [1.1]
    sorted_n = sorted(nodes)
    #return a new name of gene by selecting the last and add +1
    name = sorted_n[-1] + 1
    if name%10 == 0: name +=1
    #format to return an id
    nb = float(10**len(str(name)))
    name = formating(1. + name/nb, nb)
    return name

def formating(x, nb):
    """
    format a number
    x: number to format
    nb = cts["max_nb_of_genes"]
    """
    nb = str(int(np.log10(nb*nb)))
    s = "{0:."+ nb +"f}"
    return float(s.format(x)) #index for the node
    
def fckit(a):
    return float(str(a))
    
def verify_attributes(genome):
    for n, dic in genome.nodes(data = True):
        if dic == {}: raise AssertionError("error on the node "+str(n))
        
def verify_nbinteractions(genome, tf_lim):
    """verify that genes respect TFBS limit (tf_lim)"""
    nodes = dict(genome.nodes(data = True))
    edges = {e[0:2]:e[2] for e in genome.edges(data = True)}
    for n in nodes.keys():
        dic = nodes[n]
        if (dic["type"] == "gene" or dic["type"] == "input") and len(genome.predecessors(n)) > tf_lim:
            raise AssertionError("tf_lim is exceed for the node "+str(n))
            
def verify_protcom(genome):
    nodes = dict(genome.nodes(data = True))
    edges = {e[0:2]:e[2] for e in genome.edges(data = True)}
    for n in nodes.keys():
        dic = nodes[n]
        if dic["type"] == "protein" and len(genome.successors(n)) == 0:
            raise AssertionError("This protein has no successors: "+str(n))
        # elif dic["type"] == "complex":
            # pred = [p for p in genome.predecessors(n) if nodes[p]["type"] == "protein" and edges[(p, n)]["effect"] == "complex"]
            # if len(pred) != 2: raise AssertionError("This complex has a problem (pred: "+str(pred)+"): "+str(n))

def verify_genome(genome):
    """
    verify there is still a connection between the input and the output
    and if all nodes avec an interaction with at least one another node
    Return True if the genome is good, False if not
    """
    #interactions
    nodes = dict(genome.nodes(data = True))
    for n in nodes.keys():
        if len(genome.successors(n)) + len(genome.predecessors(n)) == 0: return False
    #inputs and outputs
    inp = [node for node, dic in genome.nodes(data = True) if dic["type"] == "input"]
    out = [node for node, dic in genome.nodes(data = True) if dic["type"] == "output"]
    for i in inp:
        if not nx.has_path(genome, i, out[0]): return False
    return True

def all_successors(genome, nodes, all_succ, init = False):
    """
    return all_successors of a node AND the node itself if init = False
    nodes = [start] where start if the parent of all successors
    """
    succ = [s for n in nodes for s in genome.successors(n)]
    new = list(set(succ) - set(nodes) - set(all_succ)) #list of new found successors
    if not init: all_succ += nodes
    if new != []: return all_successors(genome, new, all_succ, init = False)
    else: return all_succ
    
def all_succ_complexes(genome, data, nodes, all_succ, init = True):
    """return all successors of a node which are complexes"""
    succ = [s for n in nodes for s in genome.successors(n) if data[s]["type"] == "complex"]
    new = list(set(succ) - set(nodes) - set(all_succ)) #list of new found successors
    if not init: all_succ += nodes
    if new != []: return all_succ_complexes(genome, data, new, all_succ, init = False)
    else: return all_succ
    
def verify_complexes(genome, r = False):
    """return all complexes which have themself or one of their successors as a predecessor
    and verify if there are 2 predecessors for each complex"""
    nodes = dict(genome.nodes(data = True))
    errors = []
    for c in nodes:
        if nodes[c]["type"] == "complex":
            pre = genome.predecessors(c)
            suc = all_succ_complexes(genome, nodes, [c], [], init = True)
            if c in pre or c in suc or len(pre) < 2: errors += [c]
    if r: return errors
    elif errors != []: return False
    else: return True
    