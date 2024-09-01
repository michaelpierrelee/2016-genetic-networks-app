# coding: utf8

import multiprocessing
from pygenets import * 

def parameters():
    """
        For each parameter, there is a list of values which will be tested,
        such as the first value of the list of one parameter will form a set with every other first values.
        But if the user want to test few variable parameters with the others constant,
        it is not required that each list have the same length. In this case, the function
        will complete the others in function of the last element of each
        in order to that each list have the same length. For exemple:
            p = { "gen_nb": [10, 20, 30], "pop": [5, 10], "optim": [5] }
            return: p = { "gen_nb": [10, 20, 30], "pop": [5, 10, 10], "optim": [5, 5, 5] }
    """
    p = {
        "gen_nb": [1],
        "pop": [5],
        #mutations
        "optim": [45], #start of the optimization step
        "growth": [20], #end of the growth step
        "nb_mut": [0.5], #mean of nb of mutations per node per generation, nb_mut = ]0,+[
        "max_m": [10],
        "proba_optim": [
            {"m1": 0.8, "m21": 0.1, "m22": 0.1, "m31": 0, "m32": 0, "m33": 0}
        ],
        "decay": [False],
        #tournament
        "t": [0.5], #probability to keep the best, t > 0
        "score": [2], #1 = complexity score, 2 = output score
        "nb_groups": [2],
        #inputs and outputs
        "inputs": [
            {"Input": [1e-5, 0.01, 0.5, 1e3]}
        ],
        "exp_outputs": [
            [30, 30, 0, 0]
        ],
        "weight": [
            [1, 1, 1,1]
        ]
    }
    M = max([len(v) for v in p.values()]) #length longest list
    #complete the lists according to M
    for k in p.keys():
        n = len(p[k])
        if n < M: p[k] += [p[k][-1] for i in range(M - n)]
    #return
    return p, M

def new_set(i, param):
    """return a set of parameters according to the position i of each list"""
    return {
        "gen_nb": param["gen_nb"][i],
        "pop": param["pop"][i],
        #mutations
        "optim": param["optim"][i],
        "growth": param["growth"][i],
        "nb_mut": param["nb_mut"][i],
        "max_m": param["max_m"][i],
        "proba_optim": param["proba_optim"][i],
        "decay": param["decay"][i],
        #tournament
        "t": param["t"][i],
        "score": param["score"][i],
        "nb_groups": param["nb_groups"][i],
        #inputs and outputs
        "inputs": param["inputs"][i],
        "exp_outputs": param["exp_outputs"][i],
        "weight": param["weight"][i]
    }
    
def run_algo(state, param, g_init, prefix, loc):
    print "   <<< simulation " + prefix + " >>>"
    #PARAMETERS
    gen_nb = param["gen_nb"]
    pop = param["pop"]
    optim = param["optim"]
    growth = param["growth"]
    nb_mut = param["nb_mut"]
    max_m = param["max_m"]
    proba_optim = param["proba_optim"]
    decay = param["decay"]
    t = param["t"]
    score = param["score"]
    nb_groups = param["nb_groups"]
    inputs = param["inputs"]
    exp_outputs = param["exp_outputs"]
    weight = param["weight"]
    #VARIABLES
    networks = []
    pg, ct = 0, 0
    cpt = 0
    avg = [[],[]] #0 = complexity, #1 = MSE
    std = [[],[]]
    best = [[],[]] #best individual of each generation
    c = copy.deepcopy(cts())
    
    #LOOP
    for i in range(pop): new_network(networks, g_init.copy())
    for n in range(gen_nb):
        #optimization
        if n > optim: c["proba"] = proba_optim
        #mutations
        to_mutate = [g[0] for g in networks]
        networks = []
        for genome in to_mutate:
            parent = genome.copy()
            #number of mutations per generation
            m = int(nb_mut * nx.number_of_nodes(genome))
            if m == 0: m = 1
            elif m > max_m: m = max_m
            #mutations
            for i in range(m): genome = mutation(genome, c, decay)[0]
            #add the child to the list
            new_network(networks, parent)
            new_network(networks, genome)
        #incremente progress statut, displayed in the command terminal
        pg, ct = progress(cpt, gen_nb - 1, pg, ct)
        cpt += 1
        #scores
        networks = update_scores(networks, inputs, exp_outputs, weight)
        #ranking
        if n < growth: tournament(networks, score, t, nb_groups, pop)
        else: pareto(networks)
        #selection
        ranks = rank_list(networks)
        networks = selection(networks, ranks, pop, C=0)
        #monitoring
        avg[0] += [np.mean([g[1] for g in networks])]
        avg[1] += [np.mean([g[2] for g in networks])]
        std[0] += [np.std([g[1] for g in networks])]
        std[1] += [np.std([g[2] for g in networks])]
        b = sorted(networks, key=operator.itemgetter(3,2,1))[0][1:3]
        best[0] += [b[0]]
        best[1] += [b[1]]
        
    #DISPLAY
    to_log = "simulation " + prefix + "\n\n"
    #plot score evolution
    index = str(gen_nb) + "-" + str(pop) + "-" + str(cpt)
    score_evolution(loc + "scores_" + index + ".png", gen_nb, avg, std = std)
    score_evolution(loc + "best-scores_" + index + ".png", gen_nb, best)
    #plot pareto convergence
    pareto(networks, address = loc + "pareto_"+index+".png", rk = 5)
    #sorting according rank, then MSE and complexity score
    selected_net = cleanup(networks)
    selected_net = sorted(selected_net, key = operator.itemgetter(3,2,1))
    to_log += str(len(selected_net)) + " various networks generated\n\n"
    #save networks
    z = 0
    to_log += "rank\tMSE score\t\tcomplexity score\n\n"
    for i in range(len(selected_net)):
        z += 1
        genome = selected_net[i][0]
        com = round(selected_net[i][1], 3)
        out = round(selected_net[i][2],3)
        rank = selected_net[i][3]
        scores = [com, out, rank]
        to_log += str(rank) + "\t" + str(out) + "\t\t" + str(com) + "\t" + index
        if z < 5:
            index = str(rank) + "_" + str( int(round(out)) ) + "_" + str( int(round(com)) ) +"-"+str(z)
            name = loc + "curves/" + prefix + "-" + index
            if not os.path.exists(loc + "curves/"): os.makedirs(loc + "curves/")
            to_log += "\t\t(displayed)"
            #genome drawing
            #pos = layout(genome)
            pos = nx.fruchterman_reingold_layout(genome)
            pos = nx.fruchterman_reingold_layout(genome, k = 100, pos = pos, iterations = 100)
            drawing(genome, pos, name + '-structure.png', 10, 10, [-0.05,1.05], [-0.05,1.05], ax="off", scores = scores)
            #genome structure
            nx.write_yaml(genome, name + '-save.yaml')
            #dose-response curve
            i = inputs["Input"]
            inp = {"Input": np.logspace(-5, 3, 50)}
            y = test_genome(genome, inp)
            x = inp["Input"]
            ax = plt.gca()
            lab1, = plt.plot(x, y, "-r", label = "computed outputs")
            lab2, = plt.plot(i, exp_outputs, "ob", label = "expected outputs")
            ax.stem(i, exp_outputs, color="b")
            ax.set_xscale("log")
            ax.set_xlabel('Input')
            ax.set_ylabel('Output')
            plt.legend([lab1, lab2], loc = "best")
            plt.title(u"Dose–response curve for " + index)
            plt.savefig(name + "-curve.png", bbox_inches='tight')
            plt.clf()
            plt.close()
        to_log += "\n"
    #save log
    with open (loc + prefix + '-log_' + index + '.txt', 'a') as f:
        f.write(to_log)
    #return a success
    state["result"] = True
    
if __name__ == '__main__':
    param, simul = parameters()
    repeat = 1 #nb of repetitions for a simulation
    g_init = init(cts(), ["Input"])
    
    for i in range(simul):
        for j in range(repeat):
            #names
            now = datetime.datetime.now()
            prefix = "sim_" + now.strftime('%H%M-%d-%m') + "_" + str(i+1) + "_" + str(j+1)
            loc = "simulations/" + prefix + "/"
            if not os.path.exists(loc): os.makedirs(loc)
            #run
            manager = multiprocessing.Manager()
            state = manager.dict(result = False)  # shared state
            p = multiprocessing.Process(
                target = run_algo,
                args = (state, new_set(i, param), g_init, prefix, loc,)
            )
            p.start()
            p.join()
            #print in the interactive shell
            if state["result"]: msg = "Simulation " + str(i+1) + "/" + str(simul) + ", repeat " + str(j+1) + "/" + str(repeat) + " :: done"
            else:
                msg = "simulation interrupted"
                #raise AssertionError("")
            print msg
            