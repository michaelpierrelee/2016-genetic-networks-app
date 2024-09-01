# coding: utf-8

from . import *
from progress import progress

def scoring_output(genome, inputs, expected_outputs, weight):
    """
    inputs = {"input1": [values], ...} such as len([values]) == len(expected_outputs)
    weight for the expected_outputs 
    return the Mean Squared Error
    """
    MSE = 0
    for i in range(len(expected_outputs)):
        #reset input values
        dic = {k: inputs[k][i] for k in inputs.keys()}
        set_inputs(genome, dic)
        #solve the network
        measured = get_output(genome)
        #square
        MSE += ((measured - expected_outputs[i]) * weight[i]) ** 2
    return np.sqrt(MSE / len(expected_outputs))
    
def scoring_complexity(genome):
    """
    return the complexity score
    """
    complexity = 0
    edges = {e[0:2]:e[2] for e in genome.edges(data = True)}
    for n, dic in genome.nodes(data = True):
        if dic["type"] != "output" and dic["type"] != "gene":
            interact = [edges[(n, s)]["effect"] for s in genome.successors(n)
                        if edges[(n, s)]["effect"] != "transcript"
                       ]
            A = interact.count("activator")
            R = interact.count("repressor")
            C = interact.count("complex")
            complexity += 1.6**(A+R-1) * 1.25**A * 1.25**C
    return complexity
    
def network_scores(genomes, inputs, expected_outputs, weight):
    """
    genomes = [genome1, genome2, ...]
    see scoring_output for the other parameters
    return networks = [ [genome1, complex_score, output_score, rank], ... ]
    """
    MSE = []
    cpl = []
    networks = []
    for g in genomes:
        sys.stdout.write(".")
        sys.stdout.flush()
        MSE += [scoring_output(g, inputs, expected_outputs, weight)]
        cpl += [scoring_complexity(g)]
        networks += [[g, cpl[-1], MSE[-1], 1]]
    print 
    return networks
    
def reset_rank(networks):
    """reset the rank"""
    for i in range(len(networks)): networks[i][-1] = 1
    
def pareto(networks, ignored_com = False, biased = False, address = None, avg = None, sharing = False, C = 1, rk = 5, nb_sol = 5):
    """
    update the rank of the list networks following the Pareto's algorithm
    networks = [ [genome1, complex_score, output_score, rank], ... ]
    avg = means of scores
    biased = True if the output score has to be promoted from the complexity score
        it means that an element dominating another on the MSE axis has an higher rank
    """
    N = len(networks)
    
    #reset pareto rank
    reset_rank(networks)
    
    #pareto score
    for i, g1 in enumerate(networks[:N-1]):
        #pareto's algorithm
        #nets[i] = g1, nets[i+j+1] == g2
        for j, g2 in enumerate(networks[i+1:N]):
            if biased:
                if g1[1] >= g2[1] and g1[2] > g2[2]:
                    #g1 dominates g2 on Y-axis
                    networks[i][3] += 2
                elif g1[1] > g2[1] and g1[2] == g2[2]:
                    #g1 dominates g2
                    networks[i][3] += 1
                elif g1[1] <= g2[1] and g1[2] < g2[2]:
                    #g2 dominates g1 on Y-axis
                    networks[i+j+1][3] += 2
                elif g1[1] < g2[1] and g1[2] == g2[2]:
                    networks[i+j+1][3] += 1
            elif ignored_com:
                if g1[2] >= g2[2]:
                    #g1 dominates g2
                    networks[i][3] += 1
                else:
                    #g2 dominates g1
                    networks[i+j+1][3] += 1
            else:
                if (g1[1] >= g2[1] and g1[2] > g2[2]) or (g1[1] > g2[1] and g1[2] >= g2[2]):
                    #g1 dominates g2
                    networks[i][3] += 1
                elif (g1[1] <= g2[1] and g1[2] < g2[2]) or (g1[1] < g2[1] and g1[2] <= g2[2]):
                    #g2 dominates g1
                    networks[i+j+1][3] += 1
        #penalties
        if avg is not None:
            if g1[1] >= avg[0]: networks[networks.index(g1)][3] += 1
            if g1[2] >= avg[1]: networks[networks.index(g1)][3] += 1
        
    #sharing function
    if sharing:
        #radius
        for i, g in enumerate(networks):
            #initialization
            if i == 0: 
                com = [g[1], g[1]]
                out = [g[2], g[2]]
            #min and max for complex_score
            if g[1] < com[0]: com[0] = g[1]
            if g[1] > com[1]: com[1] = g[1]
            #min and max for output_score
            if g[2] < out[0]: out[0] = g[2]
            if g[2] > out[1]: out[1] = g[2]
        rs = np.sqrt((com[0] - com[1])**2 + (out[0] - out[1])**2) / (2*nb_sol)**(1/2)
        #rs = np.sqrt((com[0] - com[1])**2 + (out[0] - out[1])**2)
        #sharing
        if rs != 0:
            for i, g1 in enumerate(networks):
                s = 0
                for g2 in networks:
                    if g2 != g1:
                        shij = 1 - np.sqrt((g1[1] - g2[1])**2 + (g1[2] - g2[2])**2) / rs
                        if shij >= 0: s += shij * C
                networks[i][3] += s / N
            
    #plot
    if address != None:
        bad_MSE = [m[2] for i, m in enumerate(networks) if networks[i][3] > rk]
        bad_cpl = [m[1] for i, m in enumerate(networks) if networks[i][3] > rk]
        good_MSE = [m[2] for i, m in enumerate(networks) if networks[i][3] <= rk]
        good_cpl = [m[1] for i, m in enumerate(networks) if networks[i][3] <= rk]
        
        ax = plt.gca()
        ax.set_xscale('log', basex=2)
        ax.set_yscale('log', basey=2)
        ax.set_ylabel("Mean Squarred Error")
        ax.set_xlabel("Complexity")
        lab1, = plt.plot(bad_cpl, bad_MSE, "or", label="genomes with a rank > "+str(rk))
        lab2, = plt.plot(good_cpl, good_MSE, "ob", label="genomes with a rank <= "+str(rk))
        plt.legend([lab1, lab2,])
        plt.savefig(address, bbox_inches='tight')
        plt.clf()
        plt.close()
        
def tournament(networks, score, t, P, S):
    """
    algorithm derived from n-aire tournament
    update the rank in the list networks according only one score whose
    the position in [[genome, score1, score2, rank], ...] is given by 'score'
    t = probability to keep the best, t > 0
    P = number of groups
    S = number of individuals which have to be chosen in networks
    The list networks is updated such as the selected elements have a rank = 1
    and the others have a rank = 2
    """
    N = len(networks)
    np = N/P #group size
    sp = S/P #number of elements selected per group
    
    #reset pareto rank
    reset_rank(networks)
    
    #divide networks in sub groups which are sorted according to the specified score
    random.shuffle(networks)
    b = [sorted(networks[i : i+np], key = operator.itemgetter(score)) for i in range(0, N, np)]
    #correct the last group which can be smaller than the others
    if len(b[-2]) > len(b[-1]):
        b[-2] += b[-1]
        del b[-1]
    #select the firsts in each group with a probability t
    diff = S - len(b) * sp
    c = []
    for k, net in enumerate(b):
        #because the number of selected individuals in each group is not enough
        #to have len(c) == N/2, we will select an additionnal other element in the last groups
        if k >= N/np - diff: end = sp + 1
        else: end = sp
        #selection
        select = []
        i = 0
        # j = 0
        while len(select) != end:
            if random.uniform(0,1) < t:
                #we select the best
                select += [networks.index(net[i])]
                del net[i]
            if i >= len(net) - 1: i = 0
            else: i += 1
            # j += 1
        c += select
    #return the list networks with the rank updated:
    for i in range(N):
        if i in c: networks[i][-1] = 1
        else: networks[i][-1] = 2
        
def rank_list(networks):
    """
    return a dict sorting genomes by rank such as
        {1: [g1,g2], 2: [g3,g4], ...}
    """
    ranks = {}
    for g in networks:
        r = g[3]
        if r in ranks.keys(): ranks[r] += [g[0]]
        else: ranks[r] = [g[0]]
    return ranks
    
def update_scores(networks, inputs, expected_outputs, weight):
    """
    see network_scores() for information about parameters
    return an uptaded version of networks
    avg = means of scores
    """
    #update the scores
    to_update = [g[0] for g in networks if g[1] == 0]
    updated = [g for g in networks if g[1] != 0]
    updated += network_scores(to_update, inputs, expected_outputs, weight)
    #return
    return updated
    
def new_network(networks, genome):
    """add a new network to the list of networks"""
    networks += [[genome, 0, 0, 1]]