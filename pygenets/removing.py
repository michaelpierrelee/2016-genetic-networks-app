# coding: utf-8

from . import *
from .drawing import *

def explore_succ(genome, s, pred_c, no_link, link, nodes):
    """
    start from a node and search if its successors are linked to pred_c
    if there is a link, search in the further successors
    if there is no link, save the name of this successor in no_link and return it
    """
    succ = genome.successors(s)
    for n in succ:
        pred_succ = genome.predecessors(n)
        if pred_c in pred_succ and nodes[n]["type"] == "complex":
            link += [n]
            link, no_link = explore_succ(genome, n, pred_c, no_link, link, nodes)
        else: no_link += [[s, n]]
    return link, no_link

def remove_succ_complex(genome, c, edges, nodes, tf_lim):
    """
    remove a complex c and its successors
    /!\ need to be used only in the functions remove_
    """
    if nodes[c]["type"] == "input": raise AssertionError("Thats an input")
    succ_c = genome.successors(c)
    pred_c = genome.predecessors(c)
    pred_succ_c = [n for s in succ_c for n in genome.predecessors(s)] #predecessors of successors of the complex
    # print c
    # print "succ_c\t", succ_c
    # print "pred_c\t", pred_c
    # print "pred_succ_c\t", pred_succ_c
    done = []
    to_remove = []
    added = []
    for p in pred_c:
        if succ_c != [] and not p in pred_succ_c:
            #if the complex has a target,
            #create a link between the remaining predecessor of the complex and the target
            for s in succ_c:
                if nodes[s]["type"] == "gene" or nodes[s]["type"] == "output":
                    pred_s = [n for n in genome.predecessors(s) if n != c]
                    if p not in pred_s and len(pred_s) + 1 <= tf_lim:
                        #add a PGI between the pred (a protein) and the succ (a gene) if it is possible
                        e = edges[(c, s)]
                        add_PGI(genome, p, s, e["effect"], e["K"], e["hill"])
                        edges[(p, s)] = copy.deepcopy(e)
                        added += [(p, s)]
                    elif p not in pred_s:
                        #if there are too many proteins which already interact with the gene
                        #we remove a priori the predecessor
                        to_remove += [(p, s)]
                elif len(pred_c) == 1:
                    #get back stoechiometric coeff
                    try: pred_s = [n for n in genome.predecessors(s) if n != c][0]
                    except:
                        print c, s, p, pred_c, succ_c, pred_succ_c, done
                        print genome.predecessors(s)
                        print genome.successors(s)
                        pos = nx.fruchterman_reingold_layout(genome, k = 100, iterations = 100)
                        drawing(genome, pos, 'img/error.png', 10, 10, [-0.05, 1.05], [-0.05, 1.05])
                        raise
                    coeff1, coeff2 = nodes[s]["stoech"][c], nodes[s]["stoech"][pred_s]
                    #replace c by p
                    add_PPI(genome, p, s)
                    nx.set_node_attributes(genome, "stoech", { s: {p: coeff1, pred_s: coeff2} })
                elif not s in done:
                    done += [s]
                    #if the successor is a complex
                    #recreate the complex but only to connect predecessors to this successor
                    try: pred_s = [n for n in genome.predecessors(s) if n != c][0]
                    except:
                        print c, s, p, pred_c, succ_c, pred_succ_c, done
                        print genome.predecessors(s)
                        print genome.successors(s)
                        pos = nx.fruchterman_reingold_layout(genome, k = 100, iterations = 100)
                        drawing(genome, pos, 'img/error.png', 10, 10, [-0.05, 1.05], [-0.05, 1.05])
                        raise
                    ##recreate
                    coeff1, coeff2 = nodes[c]["stoech"][pred_c[0]], nodes[c]["stoech"][pred_c[1]]
                    kon, koff = nodes[c]["kon"], nodes[c]["koff"]
                    degrad = nodes[c]["degrad"]
                    lab = add_complex(genome, pred_c[0], pred_c[1], coeff1, coeff2, degrad, kon, koff, 1)
                    #reconfigure successor
                    add_PPI(genome, lab, s)
                    coeff1, coeff2 = nodes[s]["stoech"][c], nodes[s]["stoech"][pred_s]
                    nx.set_node_attributes(genome, "stoech", { s: {lab: coeff1, pred_s: coeff2} })
        elif succ_c != []:
            #we search the successors of successors of the complex which have no link with p
            #then we remove those with a link and we create a link between the others and pre_c
            link, no_link = explore_succ(genome, c, p, [], [], nodes)
            #print link, no_link
            for nol in no_link:
                if nodes[nol[1]]["type"] == "gene" or nodes[nol[1]]["type"] == "output":
                    pred_nol1 = genome.predecessors(nol[1])
                    if p not in pred_nol1 and len(pred_nol1) + 1 <= tf_lim:
                        #add a PGI between the pred (a protein) and the succ (a gene) if it is possible
                        e = edges[(nol[0], nol[1])]
                        add_PGI(genome, p, nol[1], e["effect"], e["K"], e["hill"])
                        edges[(p, nol[1])] = copy.deepcopy(e)
                        added += [(p, nol[1])]
                    elif p not in pred_nol1:
                        #if there are too many proteins which already interact with the gene
                        #we remove a priori the predecessor
                        to_remove += [(p, nol[1])]
                elif len(pred_c) == 1:
                    #get back stoechiometric coeff
                    try: pred_nol0 = [n for n in genome.predecessors(nol[1]) if n != nol[0]][0]
                    except:
                        print nol, c, p, pred_c
                        print genome.predecessors(nol[1])
                        print genome.successors(nol[1])
                        pos = nx.fruchterman_reingold_layout(genome, k = 100, iterations = 100)
                        drawing(genome, pos, 'img/error.png', 10, 10, [-0.05, 1.05], [-0.05, 1.05])
                        raise
                    coeff1, coeff2 = nodes[nol[1]]["stoech"][nol[0]], nodes[nol[1]]["stoech"][pred_nol0]
                    #replace c by p
                    add_PPI(genome, p, nol[1])
                    nx.set_node_attributes(genome, "stoech", { nol[1]: {p: coeff1, pred_nol0: coeff2} })
                elif not nol[1] in done:
                    done += [nol[1]]
                    #if the nol1 is a complex
                    #recreate the complex but only to connect predecessors to nol1
                    try: pred_nol0 = [n for n in genome.predecessors(nol[1]) if n != nol[0]][0]
                    except:
                        print no1, c, p, pred_c
                        print genome.predecessors(nol[1])
                        print genome.successors(nol[1])
                        pos = nx.fruchterman_reingold_layout(genome, k = 100, iterations = 100)
                        drawing(genome, pos, 'img/error.png', 10, 10, [-0.05, 1.05], [-0.05, 1.05])
                        raise
                    ##recreate
                    coeff1, coeff2 = nodes[c]["stoech"][pred_c[0]], nodes[c]["stoech"][pred_c[1]]
                    kon, koff = nodes[c]["kon"], nodes[c]["koff"]
                    degrad = nodes[c]["degrad"]
                    lab = add_complex(genome, pred_c[0], pred_c[1], coeff1, coeff2, degrad, kon, koff, 1)
                    #reconfigure successor
                    add_PPI(genome, lab, nol[1])
                    coeff1, coeff2 = nodes[nol[1]]["stoech"][nol[0]], nodes[nol[1]]["stoech"][pred_nol0]
                    nx.set_node_attributes(genome, "stoech", { nol[1]: {lab: coeff1, pred_nol0: coeff2} })
            for lin in link:
                if nodes[lin]["type"] == "complex":
                    genome.remove_node(lin)
                    done += [lin]
                    

    #removing
    if succ_c == []: to_remove += [(p, c) for p in pred_c]
    genome.remove_node(c)
    #if a successor has be linked to a predecessor but not to another
    common = [ (p1, p2, s1) for p1, s1 in to_remove for p2, s2 in added if s1 == s2]
    #randomly invert the two interactions
    for p1, p2, s in common:
        if toss():
            e = edges[(p2, s)]
            add_PGI(genome, p1, s, e["effect"], e["K"], e["hill"])
            genome.remove_edge(p2, s)
            to_remove += [(p2, s)]
            del to_remove[to_remove.index((p1, s))]
    #delete the predecessors no linked to a node (but not in the case of a pumping complex)
    for p, s in to_remove:
        if genome.has_node(p):
            if (nodes[p]["type"] == "gene" and len(genome.successors(p)) == 0) or nodes[p]["type"] == "protein":
                remove_pred(genome, p)
            elif nodes[p]["type"] == "complex" and len(genome.successors(p)) == 0:
                pred_succ_p = [v for u in genome.predecessors(p) for v in genome.successors(u)]
                if len(pred_succ_p) == 2:
                    remove_pred(genome, p)
    try:
        verify_attributes(genome)
        verify_nbinteractions(genome, tf_lim)
    except:
        print c
        pos = nx.fruchterman_reingold_layout(genome, k = 100, iterations = 100)
        drawing(genome, pos, 'img/error.png', 10, 10, [-0.05, 1.05], [-0.05, 1.05])
        raise

def remove_pred(genome, node):
    """
    remove predecessors of a node if it has only one successor (the node) without other considerations
    """
    nodes = dict(genome.nodes(data = True))
    pred_p = [node]
    purge = True
    while purge:
        purge = False #if there are no other deletions at this round, stop the loop
        #get back data for the round n+1
        pred = pred_p
        pred_p = [p for n in pred for p in genome.predecessors(n) if nodes[p]["type"] != "input"]
        #round n
        for p in pred:
            sp = [n for n in genome.successors(p) if n not in pred] #if n interacts with a node of pred
            if len(sp) == 0:
                genome.remove_node(p)
                purge = True
            #cause the n interacting with a n of pred is also in pred_p
            #or if len(sp) == 1, we don't have to remove its predecessors
            pred_p = list(set(pred_p) - set([p])) 
    
def remove_gene(genome, gene, tf_lim):
    """
    remove a gene and its protein
    replace it by an interaction between its predecessors and its successors
    such as every successors will be connected to every predecessors of the set gene/protein
    tf_lim = maximal numbers of proteins which can interact with a same gene
    """
    protein = fckit(gene + 1)
    nodes = dict(genome.nodes(data = True))
    edges = {e[0:2]:e[2] for e in genome.edges(data = True)}
    if nodes[gene]["type"] == "input" or nodes[gene]["type"] == "input": raise AssertionError("removing"+str(gene))
    out = [node for (node, dic) in genome.nodes(data = True) if dic["type"] == "output"]
    save = genome.copy()
    
    verify_attributes(save)
    verify_nbinteractions(save, tf_lim)
    
    #remove the protein and the gene
    #NB: a complex can be at the same time in succ and in pred,
    #but it will be removed by the condition "else" of the loop "for s in succ"
    succ = list(set(genome.successors(protein)) - set([gene, protein]))
    pred = list(set(genome.predecessors(gene)) - set([gene, protein]) - set(succ))
    genome.remove_node(protein)
    
    verify_attributes(genome)
    verify_nbinteractions(genome, tf_lim)
    
    #connections between the predecessors and the successors
    N = len(pred)
    to_remove = [gene]
    i = -1
    for s in succ:
        i += 1
        if genome.has_node(s):
            gs = edges[(protein, s)]
            pred_s = genome.predecessors(s)
            to_choose = list(set(pred) - set(pred_s))
            verify_nbinteractions(genome, tf_lim)
        
        if genome.has_node(s) and nodes[s]["type"] != "complex":
            #make the links between the randomly chosen predecessors and the successor "s"
            nb = tf_lim - len(pred_s)
            if nb > len(to_choose): nb = len(to_choose)
            chosen = random.sample(to_choose, nb)
            for ch in chosen:
                if ch not in pred_s:
                    add_PGI(genome, ch, s, gs["effect"], gs["K"], gs["hill"])
                    edges[(ch, s)] = copy.deepcopy(gs)
            #for the others, they are removed if there are no other succ and if they are not pumping complex
            not_chosen = list(set(to_choose) - set(chosen))
            to_remove += [ch for ch in not_chosen]
            
            try:
                verify_attributes(genome)
                verify_nbinteractions(genome, tf_lim)
            except:
                print i, s, gene, protein, pred, succ
                print azert
                print genome.predecessors(s)
                print chosen
                raise
        elif genome.has_node(s):
            #it is a complex, it will be removed and its successors will be linked to the pred of the gene
            #find the successors of the complex which are not complexes
            succ_s = genome.successors(s)
            genes = [(s,v) for v in succ_s if nodes[v]["type"] != "complex" and nodes[v]["type"] != "protein"]
            run = True
            i = 0
            while run:
                i += 1
                if i > 1000: raise AssertionError("Infinite loop")
                current = succ_s
                succ_s = [v for u in current for v in genome.successors(u) if nodes[v]["type"] == "complex"]
                genes += [(u,v) for u in current for v in genome.successors(u)
                        if nodes[v]["type"] != "complex" and nodes[v]["type"] != "protein"
                        and (u,v) not in genes
                         ]
                if succ_s == []: run = False
                
            #remove the complex
            remove_succ_complex(genome, s, {e[0:2]:e[2] for e in genome.edges(data = True)}, dict(genome.nodes(data = True)), tf_lim)
            
            try:
                verify_attributes(genome)
                verify_nbinteractions(genome, tf_lim)
            except:
                print i, s, gene, protein, pred, succ
                print genome.predecessors(s)
                raise
            #make interactions
            for c, g in genes:
                gg = edges[(c, g)]
                pred_g = genome.predecessors(g)
                to_choose = list(set(pred) - set(pred_g))
                #choose some pred and make the links
                nb = tf_lim - len(pred_g)
                if nb > len(to_choose): nb = len(to_choose)
                chosen = random.sample(to_choose, nb)
                for ch in chosen:
                    if ch not in pred_g:
                        try: add_PGI(genome, ch, g, gg["effect"], gg["K"], gg["hill"])
                        except:
                            print ch, g, gg
                            print s, genes, to_choose
                            raise
                        edges[(ch, g)] = copy.deepcopy(gg)
                        try:
                            verify_attributes(genome)
                            verify_nbinteractions(genome, tf_lim)
                        except:
                            print i, "s", s, "c", c, "p", p, "g", g
                            print "gene", gene, "protein", protein, pred, succ
                            print genome.predecessors(s)
                            raise
                #for the others, they are removed if there are no other succ and if they are not pumping complex
                not_chosen = list(set(to_choose) - set(chosen))
                to_remove += [ch for ch in not_chosen]
            try:
                verify_attributes(genome)
                verify_nbinteractions(genome, tf_lim)
            except:
                print i, s, gene, protein, pred, succ
                print genes
                pos = nx.fruchterman_reingold_layout(genome, k = 100, iterations = 100)
                drawing(genome, pos, 'img/error.png', 10, 10, [-0.05, 1.05], [-0.05, 1.05])
                raise
   
    #remove the pred
    for p in to_remove:
        if genome.has_node(p):
            if nodes[p]["type"] == "gene" and len(genome.successors(p)) == 0:
                remove_pred(genome, p)
            elif nodes[p]["type"] == "complex" and len(genome.successors(p)) == 0:
                pred_succ_p = [v for u in genome.predecessors(p) for v in genome.successors(u)]
                if len(pred_succ_p) == 2: remove_pred(genome, p)
        try:
            verify_attributes(genome)
            verify_nbinteractions(genome, tf_lim)
        except:
            print p, gene, protein, pred, succ
            print genome.predecessors(s)
            raise
            
    try: verify_protcom(genome)
    except:
        print gene, pred, succ
        pos = nx.fruchterman_reingold_layout(genome, k = 100, iterations = 100)
        drawing(genome, pos, 'img/error.png', 10, 10, [-0.05, 1.05], [-0.05, 1.05])
        raise

def remove_complex(genome, com, tf_lim):
    """
    remove a complex
    """
    nodes = dict(genome.nodes(data = True))
    edges = {e[0:2]:e[2] for e in genome.edges(data = True)}
    ##if com has successors, we link the predecessors of com to its successors
    #search neighbors of the complex
    succ = [n for n in genome.successors(com)]
    pre = [n for n in genome.predecessors(com)]
    #save interactions of the complex
    interactions = {s: edges[(com, s)] for s in succ}
    stoech = nodes[com]["stoech"]
    degrad = nodes[com]["degrad"]
    #remove the complex if it is not involved in the formation of another
    #or recreate it to only keep the connection between its 2 predecessors and the successor
    remove_succ_complex(genome, com, edges, nodes, tf_lim)
    ##if com doesn't have successors, remove parents which are only connected to the removed complex
    ##but only in the pumping complex-like configuration
    pre = [p for p in pre if genome.has_node(p)]
    pre_p = [genome.predecessors(p) for p in pre]
    suc_p = [genome.successors(p) for p in pre]
    for i, p in enumerate(pre):
        pp = pre_p[i] #if len(pp) > 1 => it is a complex
        sc = suc_p[i] #if len(sc) == 0 => it can be a pumping complex or a solitary protein
        if len(sc) == 0 and len(pp) == 1:
            remove_gene(genome, pp[0], tf_lim) #pp[0] is a prot
        elif len(sc) == 0 and len(pp) == 2: #it is a complex
            pp = [u for u in pp if genome.has_node(u)]
            suc_pp = [v for u in pp for v in genome.successors(u)]
            #if len == 2 so the complex is useless, if not it is kept because it's a pumping complex
            if len(suc_pp) == 2: remove_complex(genome, p, tf_lim)