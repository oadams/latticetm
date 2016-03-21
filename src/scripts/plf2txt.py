#!/usr/bin/python3
# -*- coding: UTF-8 -*-

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
exclude.add(u"¿")
exclude.add(u"¡")

def remove_punctuation(line):
    return "".join(c for c in line if c not in exclude)

def jointly_process_files(corpus_root, prefix, num_sents):
    """ Takes a corpus prefix and performs joint processing of (1) The
    word-lattices (2) The Spanish gold LDC transcription (3) The English
    transcription.

    (1) The LF lattices are converted to a format that can be read by latticelm.
    Sometimes there will be empty lattices. In this case, remove the
    corresponding transcription lines.

    (2) and (3) have punctuation removed and are lowercased.
    """

    plf_path = corpus_root + "/plf/" + prefix + ".es"
    gold_path = corpus_root + "/ldc/" + prefix + ".es"
    trans_path = corpus_root + "/ldc/" + prefix + ".en"

    with open(plf_path) as plf_file:
        plf_lines = plf_file.readlines()
    with open(gold_path) as gold_file:
        gold_lines = gold_file.readlines()
    with open(trans_path) as trans_file:
        trans_lines = trans_file.readlines()

    if (len(plf_lines) != len(gold_lines)) or (len(plf_lines) != len(trans_lines)):
        raise Exception("Corpora components of different lengths:\n\t" +
                        "plf_lines: %d, gold_lines: %d, trans_lines: %d" %
                            (len(plf_lines), len(gold_lines),
                            len(trans_lines)))

    lattices = [from_plf(plf) for plf in plf_lines]
    gold_lines = [remove_punctuation(line).lower() for line in gold_lines]
    trans_lines = [remove_punctuation(line).lower() for line in trans_lines]

    assert len(lattices) == len(gold_lines)
    assert len(lattices) == len(trans_lines)

    plf_path+=".n%d" % num_sents
    gold_path+=".n%d" % num_sents
    trans_path+=".n%d" % num_sents

    with open(plf_path+".lat", "w") as lat_file, open(gold_path+".nopunc.lower", "w") as gold_file, open(trans_path+".nopunc.lower", "w") as trans_file:
        for i in range(len(lattices)):
            if lattices[i] == []:
                continue
            if gold_lines[i] == "\n" or trans_lines[i] == "\n":
                continue

            if i != 0:
                print("", file=lat_file)
            for line in lattices[i]:
                print(line, file=lat_file)

            print(gold_lines[i], file=gold_file, end="")
            print(trans_lines[i], file=trans_file, end="")

            if i == num_sents-1:
                break

if __name__ == "__main__":
    corpus_root = sys.argv[1]
    prefix = sys.argv[2]
    num_sents = int(sys.argv[3])

    jointly_process_files(corpus_root, prefix, num_sents)
