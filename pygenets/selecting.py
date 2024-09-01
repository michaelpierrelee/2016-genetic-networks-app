# coding: utf8

from . import *

def selection(networks, ranks, N, C=0):
    """
    return a reduced version of the list of networks without them with a high rank
    ranks = {1: [g1,g2], 2: [g3,g4], ...}
    N = number of genomes on 1 generation. After mutation, there are a priori 2N genomes.
    networks = [ [genome1, complex_score, output_score, rank], ... ]
    C = number of elements to remove if len(networks) < 2*N
    """
    eff_N = len(networks)
    #if the number of networks is under 2N, we only keep eff_N - C networks
    if eff_N < 2*N: N = eff_N - C
    #get N first ranks
    len_r = {k: len(ranks[k]) for k in ranks.keys()}
    sort_r = sorted(len_r.keys())
    iterat = iter(sort_r)
    n = 0
    i = 0
    while n < N:
        r = next(iterat)
        n += len_r[r]
        i += 1
    #select genomes
    if n == N:
        #get selected genomes from ranks according to ranks
        to_select = [g for e in range(i) for g in ranks[sort_r[e]]]
    else:
        #get selected genomes before the last rank
        to_select = [g for e in range(i-1) for g in ranks[sort_r[e]]]
        to_choose = [g for g in ranks[r]]
        to_select += random.sample(to_choose, N - len(to_select))

    #return genomes
    return [g for g in networks if g[0] in to_select]
    
def cleanup(networks):
    """
    remove the redundancies (= two identical networks)
    networks = [ [genome1, complex_score, output_score, rank], ... ]
    """
    selected_net = []
    for i, g1 in enumerate(networks):
        nodes = sorted(g1[0].nodes(data = True))
        edges = sorted(g1[0].edges(data = True))
        other = False
        for g2 in networks[i+1:]:
            n = nodes == sorted(g2[0].nodes(data = True))
            e = edges == sorted(g2[0].edges(data = True))
            c = g1[1] == g2[1]
            o = g1[2] == g2[2]
            if n and e and c and o: other = True
        if not other: selected_net += [g1]
    return selected_net