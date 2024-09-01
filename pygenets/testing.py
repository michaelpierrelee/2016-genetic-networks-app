# coding: utf8

from . import *

def test_genome(genome, inputs):
    """
    return outputs of each input value of a genome
    inputs = {"Input1": [value11, ...], "Input2": [value21, ...]}
    """
    N = len(inputs[inputs.keys()[0]])
    test = genome.copy()
    values = []
    for i in range(N):
        #reset inputs
        inp = {k: inputs[k][i] for k in inputs.keys()}
        set_inputs(test, inp)
        #test dose-response effect
        values += [get_output(test)]
    return values