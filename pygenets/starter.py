# coding: utf-8

from . import *

def cts():
    """
    kon and koff: Corzo (2006) Time, the Forgotten Dimension of Ligand Binding Teaching.
    """
    cts = {
    "input_qty": { "default":0 },
    "output": { "degradation": 0.01, "age": 1},
    "pmax": { "default": 0.1, "domain": [1e-7, 1e-1], "variation": "log" },
    "pmin": { "default": 0.001, "domain": [1e-9, 1e-3], "variation": "log" },
    "K": { "default": 0.01, "domain": [1e-4,1e3], "variation": "log" },
    "hill": { "default": 1, "domain": [1,4], "variation": 0.5 },
    "stoech": { "default": 1, "domain": [1,5], "variation": 1 },
    "degrad": { "default": 0.01, "domain": [1e-3,1e-1], "variation": "log" },
    "max_nb_of_genes": 10,
    "TFBS_limit": 3,
    "proba":{ #the sum has to be 1 and m1>m2>m31>m32>m33
            "m1": 0.5, #modify a parameter
            "m21": 0.20, #remove a protein-gene interaction
            "m22": 0.1, #remove a gene or a complex
            "m31": 0.1, #add a new gene
            "m32": 0.1, #create a new protein-gene interaction
            "m33": 0 #create a new protein-protein interaction
        },
    "kon": { "default": 1e4, "domain": [1e4, 1e8], "variation": "log" },
    "koff": { "default": 1e5, "domain": [1e-6, 1e5], "variation": "log" },
    }
    
    return cts

def init(cts, inputs):
    """
    define a basic graph
    inputs = list of inputs
    """
    a = formating(1,10)
    genome = nx.DiGraph()
    #add intermediate nodes
    add_gene(genome, 1.1, cts["pmax"]["default"], cts["pmin"]["default"], cts["degrad"]["default"], 1)
    #add inputs
    for i in inputs:
        genome.add_node(i)
        nx.set_node_attributes(genome, "type", {
                i: "input",
            })
        nx.set_node_attributes(genome, "input_qty", {
                i: cts["input_qty"]["default"]
            })
        nx.set_node_attributes(genome, "age", {
                i: 1
            })
        add_PGI(genome, i, 1.1, "activator", cts["K"]["default"], cts["hill"]["default"])
    #add output
    genome.add_node("Output")
    nx.set_node_attributes(genome, "type", {
            "Output": "output"
        })
    nx.set_node_attributes(genome, "pmax", {
            "Output": cts["pmax"]["default"]
        })
    nx.set_node_attributes(genome, "pmin", {
            "Output": cts["pmin"]["default"]
        })
    nx.set_node_attributes(genome, "output_degradation", {
            "Output": cts["output"]["degradation"]
        })
    nx.set_node_attributes(genome, "age", {
            "Output": cts["output"]["age"]
        })
    add_PGI(genome, 2.1, "Output", "activator", cts["K"]["default"], cts["hill"]["default"])
    
    return genome
    
