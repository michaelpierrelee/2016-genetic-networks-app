# coding: utf-8

from . import *

def mut(variation, domain):
    """
    return a new value according to the type and the domain of variation,
    for a specified parameter in the variable "cts"
        variation = cts[parameter]["variation"]
        domain = cts[parameter]["domain"]
    """
    if variation == 1:
        return random.randint(domain[0], domain[1])
    elif variation == 0.5:
        return random.randint(domain[0]*2, domain[1]*2)/2.
    elif variation == "log":
        a, b = int(np.log10(domain[0])), int(np.log10(domain[1]))
        o = random.randint(a, b)
        u = random.randint(1, 9)
        return u * 10**o
        
def proba_mut(nb, gamma, decay):
    """
    nb = number of elements per group
    return the probability to obtain a group according to the exponential (decay) law
    such as the first group has the biggest probability and the last group the lowest in the case of the decay law
    decay = True => follow the exp decay law, if decay = False = > the exp law
    if decay = None => Uniform law <=> gamma = 0.1
    The probability of a group is also dependent on the number of elements in a group
    in order to increase the probability of groups with a bigger number of elements
    if gamma == 0, the probability is uniform for each group
    """
    N = len(nb)
    m = float(sum(nb))
    x = np.arange(1, N + 1)
    if decay is None or gamma == 0: #uniform law
        y = np.array([nb[i] / m for i in range(N)])
    elif decay: #exponential decay law
        y = np.array([(nb[i] / m) * gamma * np.exp(- x[i] * gamma) for i in range(N)]) 
    else: #exponential law
        y = np.array([(nb[i] / m) * gamma * np.exp(x[i] * gamma) for i in range(N)]) 
    y = y / sum(y)
    return list(y)
        
def choose_node(nodes, decay, P = 5, gamma = 1):
    """
    chose a node in the list according to its age
    nodes = [[node1, node1 age], ...]
    for the parameter decay, see proba_mut
    P = number of groups of ages. It means if there are 6 ages (1,2,3,4,5,6) for P=2,
    so the first group will avec the first ages (1,2,3) and the other the lasts (4,5,6)
    gamma = constant of the exponential law, cf function proba_mut
    """
    #nodes = [[n, dic["age"]] for n, dic in genome.nodes(data = True) if dic["type"] is in types]
    nodes = sorted(nodes, key = operator.itemgetter(1))
    ages = [n[1] for n in nodes]
    #find the number of different ages per group of ages
    qty = Counter(ages) 
    qty = sorted([(k, nb) for k, nb in dict(qty).iteritems()], key = operator.itemgetter(0))
    nbByAges = [n[1] for n in qty]
    #normalization of ages to produce P groups with several ages
    ages = [n[0] for n in qty] #list of ages without redundancies
    N = len(ages)
    if N < P: P = N
    np = N/P
    diff = N - P * np
    groups = [ages[i:i+np] for i in range(0, np*(P - diff), np)]
    nb = [sum(nbByAges[i:i+np]) for i in range(0, np*(P - diff), np)]
    #create the last groups with one more individual per group
    if diff > 0:
        groups += [ages[i:i+np+1] for i in range(np*(P - diff), N, np + 1)]
        nb += [sum(nbByAges[i:i+np]) for i in range(np*(P - diff), N, np + 1)]
    #give a probability to each group
    proba = proba_mut(nb, gamma, decay)
    #select a group of ages
    r = random.random()
    j = 0
    selected = 0
    a = None
    while j < P + 1:
        if r > 1 - sum(proba[:j + 1]):
            selected = groups[j]
            a = j
            j = P + 1
        else: j += 1
    #select an age in the group
    age = random.choice(selected)
    #select a node with this age
    to_choice = [n for n, a in nodes if a == age]
    return random.choice(to_choice)
    
def mutate_param(genome, cts, decay):
    """
    chose a node or an edge and modify one of its parameter (pmin, pmax or degrad for nodes, K, hill or effect for edges)
    Cannot modify the output parameters
        cts = contain constants
    """
    edges = {(u, v):dic for (u, v, dic) in genome.edges(data = True)
             if dic["effect"] == "activator"
             or dic["effect"] == "repressor"
            }
    #selection
    #there is a bigger probability to modify edges than nodes, because there is a lot of it.
    a = nx.number_of_nodes(genome)
    b = len(edges)
    if proba_uni(1, a+b, a):
        nodes = {node: dic for (node, dic) in genome.nodes(data = True) if dic["type"] != "input"}
        where_choose = [[n, dic["age"]] for n, dic in genome.nodes(data = True)
                        if dic["type"] in ["output", "complex", "protein", "gene"]
                       ]
        chosen = choose_node(where_choose, decay)
        if nodes[chosen]["type"] == "gene"or nodes[chosen]["type"] == "output":
            #select a gene and modify its pmax or min
            if toss(): mutate = "pmax"
            else: mutate = "pmin"
            m = mut(cts[mutate]["variation"], cts[mutate]["domain"])
            #verify that pmax >= pmin
            if mutate == "pmax" and m < nodes[chosen]["pmin"]:
                nx.set_node_attributes(genome, "pmin", {chosen: m})
                nx.set_node_attributes(genome, mutate, {chosen: nodes[chosen]["pmin"]})
            elif mutate == "pmin" and m > nodes[chosen]["pmax"]:
                nx.set_node_attributes(genome, "pmax", {chosen: m})
                nx.set_node_attributes(genome, mutate, {chosen: nodes[chosen]["pmax"]})
            else:
                nx.set_node_attributes(genome, mutate, {chosen: m})
        elif nodes[chosen]["type"] == "protein":
            #select a protein and modify its degradation coefficient
            m = mut(cts["degrad"]["variation"], cts["degrad"]["domain"])
            nx.set_node_attributes(genome, "degrad", {chosen: m})
        else:
            #select a complex and modify one stoechiometric coefficient or its degradation
            a = random.randint(1,4)
            if a == 1:
                mutate = "stoech"
                prot = nodes[chosen][mutate]
                if toss(): prot[prot.keys()[0]] = mut(cts[mutate]["variation"], cts[mutate]["domain"])
                else: prot[prot.keys()[1]] = mut(cts[mutate]["variation"], cts[mutate]["domain"])
            elif a == 2:
                mutate = "kon"
                prot = mut(cts[mutate]["variation"], cts[mutate]["domain"])
            elif a == 3:
                mutate = "koff"
                prot = mut(cts[mutate]["variation"], cts[mutate]["domain"])
            else:
                mutate = "degrad"
                prot = mut(cts[mutate]["variation"], cts[mutate]["domain"])
            nx.set_node_attributes(genome, mutate, {chosen: prot})
    else:
        #select an edge of a chosen node
        types = ["input", "protein", "complex"]
        where_choose = [
            [n, dic["age"]]
            for n, dic in genome.nodes(data = True)
            for s in genome.successors(n)
            if dic["type"] in types and (n, s) in edges
            ]
        node_chosen = choose_node(where_choose, decay)
        edges = [(e[0], e[1]) for e in genome.edges() if e[0] == node_chosen]
        #modify K, effect or hill parameter
        if edges != []:
            chosen = random.choice(edges)
            r = random.randint(1,3)
            if r == 1:
                nx.set_edge_attributes(genome, "K", {
                        (chosen[0], chosen[1]): mut(cts["K"]["variation"], cts["K"]["domain"])
                    })
            elif r == 2:
                nx.set_edge_attributes(genome, "hill", {
                        (chosen[0], chosen[1]): mut(cts["hill"]["variation"], cts["hill"]["domain"])
                    })
            else:
                nx.set_edge_attributes(genome, "effect", {
                        (chosen[0], chosen[1]): random.choice(["activator", "repressor"])
                    })

def mutate_remove(genome, tf_lim, decay):
    """
    remove a gene or a complex
    """
    #select a node to remove
    nodes = {node: dic for (node, dic) in genome.nodes(data = True)
             if dic["type"] == "gene"
             or dic["type"] == "complex"
            }
    if nodes == {}: return False
    types = ["gene", "complex"]
    chosen = choose_node([[n, dic["age"]] for n, dic in genome.nodes(data = True) if dic["type"] in types], decay)
    #remove
    if nodes[chosen]["type"] == "gene": remove_gene(genome, chosen, tf_lim)
    else: remove_complex(genome, chosen, tf_lim)
    
def mutate_removeInteraction(genome, tf_lim, decay):
    """
    Remove randomly selected protein-gene and protein-protein interactions
    cannot remove an interaction protein->gene if:
        - it is the only interaction for the protein, even if there is a feedback 
        - the gene is the output which has only one interaction
    """
    nodes = dict(genome.nodes(data = True))
    edges = {e[0:2]:e[2] for e in genome.edges(data = True)
             if e[2]["effect"] == "complex"
             or (nodes[e[1]]["type"] != "output"
                 and (e[2]["effect"] == "activator" or e[2]["effect"] == "repressor")
                 and len( set(genome.successors(e[0])) - set(genome.predecessors(e[0])) ) > 1
                )
             or (nodes[e[1]]["type"] == "output"
                 and len(genome.predecessors(e[1])) > 1
                 and len( set(genome.successors(e[0])) - set(genome.predecessors(e[0])) ) > 1
                )
            }
            
    if len(edges) >= 1:
        #selection of a node
        where_choose = [e[0] for e in edges] #select all first possible nodes
        where_choose = list(set(where_choose) & set(where_choose)) #remove redundancies
        to_choose = [[e, nodes[e]["age"]] for e in where_choose] #prepare the list
        chosen = choose_node(to_choose, decay) #choose a node
        #selection of an edge
        target = random.choice([e for e in edges if e[0] == chosen])
        #test the operation on a copy
        if edges[target]["effect"] == "complex":
            remove_complex(genome, target[1], tf_lim)
        else:
            genome.remove_edge(target[0], target[1])

def mutate_addinteract(genome, cts, decay, protein = None, com = True):
    """
    Give a function to a protein. Select a target node and assign an interaction.
    If it is a gene, a PGI is added.
    If it is a protein, a complex is added (see mutate_addcomplex()).
    If com = False, so the available targets are only genes.
    But it can have a case where all genes and outputs are not available anymore (when TFBS limit is reached).
    So, we start the remove_complex() function on this added complex,
    but not in the case where it is a pumping complex.
    protein = node label of the protein where starts the interaction
    Return False if the function failed (so if we start a removing function)
        protein: specified node label which is the starting point
        com: True if "protein" can form a new complex with another, False if it can be only
        a transcription factor
    """
    nodes = {}
    proteins = []
    inputs = []
    genes = []
    for (node, dic) in genome.nodes(data = True):
        nodes[node] = dic
        if dic["type"] == "protein" or dic["type"] == "complex":
            proteins += [[node, dic["age"]]]
        elif dic["type"] == "input":
            inputs += [[node, dic["age"]]]
        elif len(genome.predecessors(node)) < cts["TFBS_limit"]:
            genes += [[node, dic["age"]]]

    #protein/input selection
    adding = True
    if protein == None:
        if not com: protein = choose_node(proteins + inputs, decay)
        else: protein = choose_node(proteins, decay)
        adding = False #because mutate_addinteract() is not used by another mutate_ function

    #cleaning of lists: we cannot target a successor of the protein
    succ = [[n, nodes[n]["age"]] for n in genome.successors(protein)]
    proteins = [p for p in proteins if p != [protein, nodes[protein]["age"]] and p not in succ]
    if adding:
        #if the protein has been added by a mutate_ function,
        #so it must not make an interaction with its own gene,
        #if not an independent network would be created
        pred = [[n, nodes[n]["age"]] for n in genome.predecessors(protein)]
        genes = [g for g in genes if g not in succ and g not in pred]
    else:
        genes = [g for g in genes if g not in succ]

    #target selection
    if com and genes+proteins != []: target = choose_node(genes + proteins, decay)
    elif genes != []: target = choose_node(genes, decay)
    else: target = None

    #add interaction
    if target != None:
        if nodes[target]["type"] == "gene" or nodes[target]["type"] == "output":
            #it will be a transcription factor
            effect = random.choice(["activator", "repressor"])
            K = mut(cts["K"]["variation"], cts["K"]["domain"])
            hill = mut(cts["hill"]["variation"], cts["hill"]["domain"])
            add_PGI(genome, protein, target, effect, K, hill)
        elif com:
            #it will be a complex
            return mutate_addcomplex(genome, cts, decay, protein1 = protein, protein2 = target)
    elif adding:
        #it has to be a transcription factor but there is no available genes
        #so we remove the gene or the complex
        if nodes[protein]["type"] == "protein":
            remove_gene(genome, fckit(protein - 1.), cts["TFBS_limit"])
            return False
        elif nodes[protein]["type"] == "complex":
            pre = genome.predecessors(protein)
            succ_pre = [s for p in pre for s in genome.successors(p)]
            #if len = 2, new complex is useless, remove it. If not, it is a pumping complex
            if len(succ_pre) == 2:
                remove_complex(genome, protein, cts["TFBS_limit"])
                return False
    
def mutate_addgene(genome, cts, decay):
    """
    Add a new gene and its protein. A function will be given to it.
    See mutate_addinteract() to have more information about this last case.
    Return False if the function failed
        cts = contain constants
    """
    gene = gene_name(genome)
    protein = fckit(1. + gene)
    #add the new gene
    pmax = mut(cts["pmax"]["variation"], cts["pmax"]["domain"])
    pmin = mut(cts["pmin"]["variation"], cts["pmin"]["domain"])
    degrad = mut(cts["degrad"]["variation"], cts["degrad"]["domain"])
    add_gene(genome, gene, pmax, pmin, degrad, 1)
    #give a function to the new protein
    mutate_addinteract(genome, cts, decay, protein = protein)

def mutate_addcomplex(genome, cts, decay, protein1 = None, protein2 = None):
    """
    Add a complex from defined proteins or not. A function is given to this complex.
    It can be only a transcription factor, to avoid the formation of an infinite complexification.
    See mutate_addinteract() to have more information 
    Return False if the function failed
        cts = contain constants
    """
    #choice of proteins if they are not defined
    nodes = dict(genome.nodes(data = True))
    complexes = [node for node, dic in genome.nodes(data = True)
                 if dic["type"] == "complex" or dic["type"] == "protein"
                ]
    if complexes == []: return False
    if protein1 == None or protein1 == protein2: #chose a first protein
        c = list(set(complexes) - set([protein2]))
        to_choose = [[n, nodes[n]["age"]] for n in c]
        protein1 = choose_node(to_choose, decay)
    if protein2 == None or protein1 == protein2: #chose a second protein
        #we cannot chose 2 proteins already giving a complex
        pre = [p for s in genome.successors(protein1)
               for p in genome.predecessors(s)
               if nodes[s]["type"] == "complex"
              ]
        c = list(set(complexes) - set([protein1]) - set(pre))
        to_choose = [[n, nodes[n]["age"]] for n in c]
        protein2 = choose_node(to_choose, decay)
    #choice of stoechiometric coefficients
    a = mut(cts["stoech"]["variation"], cts["stoech"]["domain"])
    b = mut(cts["stoech"]["variation"], cts["stoech"]["domain"])
    kon = mut(cts["kon"]["variation"], cts["kon"]["domain"])
    koff = mut(cts["koff"]["variation"], cts["koff"]["domain"])
    degrad = mut(cts["degrad"]["variation"], cts["degrad"]["domain"])
    #add the complex
    new_complex = add_complex(genome, protein1, protein2, a, b, degrad, kon, koff, 1)
    #give a function to this complex, it will only be transcripton factor
    #if it is in a pumping complex-like configuration, chose randomly if it will have a fct or not
    pred_succ = [v for u in genome.predecessors(new_complex) for v in genome.successors(u)]
    if len(pred_succ) <= 2 or toss():
        #if it is not in a complex-like configuration or if the toss succeeds:
        mutate_addinteract(genome, cts, decay, protein = new_complex, com = False)

def mutation(genome, cts, decay, i = 1):
    if not verify_genome(genome) or not verify_complexes(genome): raise AssertionError("yolo")
    p = cts["proba"]
    r = random.randint(0,100)/100.
    gene_nb = len([n for n, dic in genome.nodes(data = True) if dic["type"] == "gene"])
    save = genome.copy()
    m = ""
    if i > 4 or (r >= 1 - p["m1"] or (r >= (1 - p["m1"] - p["m21"] - p["m22"] - p["m31"]) and gene_nb > cts["max_nb_of_genes"])):
        #m1, modify a parameter
        mutate_param(save, cts, decay)
        m = "modify a parameter"

    elif r >= 1 - p["m1"] - p["m21"]:
        #m21, remove a protein-gene interaction
        mutate_removeInteraction(save, cts["TFBS_limit"], decay)
        m = "remove a protein-gene interaction"
        
    elif r >= 1 - p["m1"] - p["m21"] - p["m22"]:
        #m22, remove a gene or a complex
        mutate_remove(save, cts["TFBS_limit"], decay)
        m = "remove a gene or a complex"

    elif r >= (1 - p["m1"] - p["m21"] - p["m22"] - p["m31"]) and gene_nb <= cts["max_nb_of_genes"]:
        #m31, add a new gene
        #if the network is already full of genes, the default action is m1
        mutate_addgene(save, cts, decay)
        m = "add a new gene"

    elif r >= (1 - p["m1"] - p["m21"] - p["m22"] - p["m31"]- p["m32"]):
        #m32, create a new gene-protein interaction
        mutate_addinteract(save, cts, decay, com = False)
        m = "create a new gene-protein interaction"

    elif r >= (1 - p["m1"] - p["m21"] - p["m22"] - p["m31"]- p["m32"] - p["m33"]):
        #m33, create a new protein-protein interaction
        mutate_addinteract(save, cts, decay, com = True)
        m = "create a new protein-protein interaction"
    
    if not verify_genome(save) or not verify_complexes(save):
        return mutation(genome, cts, decay, i = i+1)
   
    return save, m, i