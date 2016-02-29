#!/usr/bin/python

from __future__ import print_function
import sys
import io
import string

def multiple_end_states(plf):
    """ Verifies that a PLF has only one final state."""

    end_i = len(plf)
    i = 0
    for out_edges in plf:
        for edge in out_edges:
            if i + edge[2] > end_i:
                return True
    return False

def from_plf(plf):
    """ Takes a string representing a single lattice in Python Lattice Format
    (PLF) and outputs a textual representation for use with latticelm."""

    # Secure your lattices, I guess?
    try:
        plf = eval(plf)
    except:
        return []

    assert(not multiple_end_states(plf))

    num_nodes = len(plf)+1 # Add one for the starting node.

    i = 0
    output_lattice = []
    for outgoing_edges in plf:
        for edge in outgoing_edges:
            # Converting logs to negative logs with the last argument.
            output_lattice.append("%d\t%d\t%s\t%s\t%f" % (
                    i, i+edge[2], edge[0], edge[0], edge[1]*-1))
        i += 1
    return output_lattice

exclude = set(string.punctuation)
def process_english(line):
    return "".join(c.lower() for c in line if c not in exclude)

def convert_file(es_in_fn, es_out_fn, en_in_fn, en_out_fn):
    """ Takes input and output filenames as strings. Reading PLF lattices from
    the input file outputs lattices in a format to be read by latticelm.
    Sometimes there will be empty lattices. In this case, remove the
    corresponding English line."""
    with open(es_in_fn) as input_file:
        es_input_lines = input_file.readlines()
    with io.open(en_in_fn, encoding="utf-8") as input_file:
        en_input_lines = input_file.readlines()
    with open(es_out_fn, "w") as es_out_file, io.open(en_out_fn, "w", encoding="utf-8") as en_out_file:
        i = 0
        for i in range(len(es_input_lines)):
            plf = es_input_lines[i]
            output_lattice = from_plf(plf)
            #print(output_lattice)
            #raw_input()
            if output_lattice == []:
                continue
            for line in output_lattice:
                print(line, file=es_out_file)
            print("",file=es_out_file)
            print(process_english(en_input_lines[i]), file=en_out_file, end=u"")
            #print(process_english(en_input_lines[i]))

#plf = "((('einen',1.0,1),),(('wettbewerbsbedingten',0.5,2),('wettbewerbs',0.25,1),('wettbewerb',0.25, 1),),(('bedingten',1.0,1),),(('preissturz',0.5,2),('preis',0.5,1),),(('sturz',1.0,1),),)"
#from_plf(plf)

es_in_fn, es_out_fn, en_in_fn, en_out_fn = sys.argv[1:]
convert_file(es_in_fn, es_out_fn, en_in_fn, en_out_fn)
