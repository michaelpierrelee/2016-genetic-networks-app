# coding: utf-8

from . import *

X = []
Y = []

def initPlot(y):
    """prepare the plot to display the solutions proposed by fsolve"""
    global X
    global Y
    if X != []: X += [X[-1] + 1]
    else: X = [0]
    Y += [y]
    
def printPlot(address):
    """display the plot with in X-axis the iteration and in Y-axis the
    solution proposed by fsolve"""
    global X
    global Y
    ax = plt.gca()
    #ax.set_yscale('log')
    plt.plot(X, Y, "or")
    plt.savefig(address, bbox_inches='tight')
    plt.clf()
    plt.close()

def hill_eqt(proteins, t, *args):
    """
    proteins = (X, Y, ...): tuple of solutions for each concentration
    factors = {X: [X factors], Y: [Y factors], ...}: parameters of
        the differential eqt dX/dt, dY/dt... such as
        for type = "protein":
            [factors of X] =
            [
                type, pmin, pmax, degrad,
                [[TF parameters of A], [TF parameters of B], ...],
                [[complex parameters of C], [complex parameters of D], ...]
            ]
            - where "TF parameters A" are the parameters for the TF A acting on the X gene
            such as [TF parameters A] = [name, effect, K, hill]
            - and "complex parameters of C" are the parameters of the complex C
            involving the protein X and another in its formation
            such as [param] = [name A, name C, X stoech coeff, A stoech coeff, kon, koff]
        for type = "complex":
            [N factors] = [type, degrad, {stoech A}, {stoech B}, kon, koff]
            where {stoech A} = {"name":"A", "coeff":1} for example
    names = ["X", "Y", ...]
        where the position 0 corresponds to the element at the position 0 in the list "proteins"
    qty = {X:, Y:, A:, B:, ...}
        where each entry is the quantity of the element specified by the key
    """
    eqt = []
    factors, names, qty = args
    #print proteins
    #initPlot(qty["Output"])

    #update the list of quantities
    for i, p in enumerate(proteins):
        qty[names[i]] = float(p)
        # if (abs(p) > 1e100 or abs(p) < 1e-100) and p != 0:
            # raise AssertionError("too large:\n" + str(i) + "\n\n" + str(p) + "\n\n" + str(names) + "\n\n" + str(proteins)  + "\n\n" + str(qty)+ "\n\n" + str(factors))
            
    #solve the equations
    for i, p in enumerate(proteins):
        prot = factors[names[i]]
        if prot[0] == "protein":
            #dX/dt = hill - degrad - complex
            pmin = float(prot[1])
            pmax = float(prot[2])
            degrad = float(prot[3])
            act = 0.
            rep = 0.
            on = 0.
            off = 0.
            for f in prot[4]:
                q = float( qty[f[0]] )
                if q < 0: q = 0
                if f[1] == "activator":
                    act += ( q / f[2] ) ** f[3]
                else: 
                    rep += ( q / f[2] ) ** f[3]
            for f in prot[5]:
                q = float( qty[f[0]] )
                c = float( qty[f[1]] )
                print " ", f[0], names[i], f[1], 
                on += kon * p**f[2] * q**f[3]
                off += koff * c
            expr = pmin + (pmax - pmin) * (act / (1 + act)) * (1 / (1 + rep)) - degrad * p - on + off
        elif prot[0] == "complex":
            #dC/dt = kon * A * B - koff * complex - degrad
            degrad = float( prot[1] )
            kon, koff = float( prot[4] ), float( prot[5] )
            A = float( qty[prot[2]["name"]] ) ** prot[2]["coeff"]
            B = float( qty[prot[3]["name"]] ) ** prot[3]["coeff"]
            expr = kon * A * B - koff * p - degrad * p

        eqt += [expr]
            
    return tuple(eqt)

def hill_eqt2(t, proteins, a, b, c):
    args = a, b, c
    tup = hill_eqt(proteins, t, *args)
    return list(tup)

def prepare(genome):
    """
    
    """
    nodes = {node: dic for node, dic in genome.nodes(data = True)}
    prots = [n for n in nodes.keys() if nodes[n]["type"] != "gene"]
    edges = {(u,v): dic for u,v,dic in genome.edges(data = True)}
    
    proteins = []
    names = []
    factors = {}
    qty = {}
    
    for node in prots:
        n = nodes[node]
        if n["type"] == "protein" or n["type"] == "complex":
            qty[node] = 0
            proteins += [0]
            names += [node]
            if n["type"] == "protein":
                #node parameters
                gene = genome.predecessors(node)[0]
                pmin, pmax = nodes[gene]["pmin"], nodes[gene]["pmax"]
                degrad = nodes[node]["degrad"]
                #TF parameters
                pred = genome.predecessors(gene)
                TF_param = []
                for p in pred:
                    edge = edges[(p, gene)]
                    TF = [p, edge["effect"], edge["K"], edge["hill"]]
                    TF_param += [TF]
                #complex parameters
                succ = genome.successors(node)
                com_param = []
                for s in succ:
                    if nodes[s] == "complex":
                        stoech = nodes[s]["stoech"]
                        kon, koff = nodes[s]["kon"], nodes[s]["koff"]
                        com_param += [[s, stoech[node], stoech[s], kon, koff]]
                #set
                factors[node] = ["protein", pmin, pmax, degrad, TF_param, com_param]
            else:
                #parameters
                kon, koff = nodes[node]["kon"], nodes[node]["koff"]
                degrad = nodes[node]["degrad"]
                #stoechiometric coefficients
                s = nodes[node]["stoech"].keys()
                stoech_A = {}
                stoech_A["name"] = s[0]
                stoech_A["coeff"] = nodes[node]["stoech"][s[0]]
                stoech_B = {}
                stoech_B["name"] = s[1]
                stoech_B["coeff"] = nodes[node]["stoech"][s[1]]
                #set
                factors[node] = ["complex", degrad, stoech_A, stoech_B, kon, koff]
        elif n["type"] == "input":
            qty[node] = nodes[node]["input_qty"]
        else: #output
            qty[node] = 0
            proteins += [0]
            names += [node]
            #factors
            pmin, pmax = nodes[node]["pmin"], nodes[node]["pmax"]
            degrad = nodes[node]["output_degradation"]
            #TF
            pred = genome.predecessors(node)
            TF_param = []
            for p in pred:
                edge = edges[(p, node)]
                TF = [p, edge["effect"], edge["K"], edge["hill"]]
                TF_param += [TF]
            #set
            factors[node] = ["protein", pmin, pmax, degrad, TF_param, []]
    
    
    return proteins, factors, names, qty
    
def handle_warning(message, category, filename, lineno, file=None, line=None):
    pass
    
def solver(genome):
    """
    see hill_eqt() for information about parameters
    Solve the Hill's equations at the steady state
    """
    global X
    global Y
    X = []
    Y = []
    proteins, factors, names, qty = prepare(genome)
    with warnings.catch_warnings():
        warnings.showwarning = handle_warning #catch the warnings from fsolve
        f = scipy.optimize.fsolve(hill_eqt, proteins, (factors, names, qty))

    return names, f
    
def integrater(genome, address = None, a = 10, b = 10, E = 100, epsilon = 0.01):
    """
    a = end time for a time interval
    b = number of points taken in the interval [0, a]
    E = maximal number of rounds, the final end time is so equal to a * E
    epsilon = difference between the last value of an interval and the mean value on it
    (an interval is a list generated by odeint of same length as np.linspace(0, a, b))
    see hill_eqt() for information about parameters
    make the integration of Hill's equations
    Return the final values of each concentration when the difference between
    the current value and the mean on the "b" last values is inferior to epsilon
    or after E rounds (the returned final value is so equal to the mean on the last values)
    """
    proteins, factors, names, qty = prepare(genome)
    N = len(proteins)
    t = np.linspace(0, a, b)
    # if address is not None: results = np.array([proteins])
    
    results = np.array([proteins])

    e = 0
    end = False
    while not end and e < E:
        e += 1
        #integration
        with warnings.catch_warnings(): #there are errors but we don't display them
            warnings.filterwarnings('ignore')
            f, msg = scipy.integrate.odeint(hill_eqt, proteins, t, (factors, names, qty), full_output = 1)
        # if address is not None: results = np.concatenate((results,f[:-1]))
        results = np.concatenate((results,f[:-1]))
        #difference 
        end = True
        for i in range(N):
            if abs(f[-1, i] - np.mean(f[:, i])) > epsilon: end = False
        #prepare next round
        proteins = [p for p in f[-1]]

    if address is not None:
        lab = []
        plt.figure(figsize=(10*1.12, 10*1.19), dpi=100)
        for i in range(N):
            l, = plt.plot(results[:,i], label=names[i], marker=".")
            lab += [l]
        ax = plt.gca()
        ax.set_yscale("log")
        plt.legend(lab, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(address, bbox_inches='tight', dpi = 100)
        plt.clf()
        plt.close()
        return names, results
    else: return names, f
    
def integrater2(genome, address = None, a = 10, b = 10, E = 1000, epsilon = 0.001):
    """
    a = nb of last values taken to compute the mean which will be compared to the last value
    and if the difference between them is less than epsilon, so this last value will be returned
    a/b = time step
    E = maximal number of integration
    """
    proteins, factors, names, qty = prepare(genome)
    N = len(proteins)
    dt = a/float(b)
    ys = np.array([proteins])
    end = False
    
    od = ode(hill_eqt2)
    # od.set_integrator('cvode', nsteps = 500, method = 'bdf')
    od.set_integrator('dopri5', verbosity = 0) #there are errors but we don't display them
    od.set_initial_value(proteins, 0).set_f_params(factors, names, qty)
    
    while od.successful() and od.t < E*dt and not end:
        # with warnings.catch_warnings():
            # warnings.filterwarnings('ignore')
            # warnings.showwarning = handle_warning
        od.integrate(od.t + dt)
        ys = np.concatenate((ys, [od.y]))
        end = True
        for i in range(N):
            if od.t < a * dt or abs(ys[-1, i] - np.mean(ys[-a:, i])) > epsilon: end = False

    if address is not None:
        lab = []
        plt.figure(figsize=(10*1.12, 10*1.19), dpi=100)
        for i in range(N):
            l, = plt.plot(ys[:,i], label=names[i], marker=".")
            lab += [l]
        ax = plt.gca()
        ax.set_yscale("log")
        plt.legend(lab, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(address, bbox_inches='tight', dpi = 100)
        plt.clf()
        plt.close()
            
    return names, ys[-a:]
    
def get_output(genome):
    """
    return the output value of a genome
    It is the unique function of this file which has to keep its name
    """
    names, solutions = integrater(genome)
    solutions = [np.mean(solutions[:, i]) for i in range(len(names))]
    output = 0
    for i, n in enumerate(names):
        if n == "Output": output = solutions[i]
    return output