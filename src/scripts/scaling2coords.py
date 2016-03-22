#!/usr/bin/python3

import sys

results_path = sys.argv[1]

with open(results_path) as results_f:
    lines = results_f.readlines()

for line in lines:
    split = line.split()
    print("(%s,%s)" % (split[0],split[1]))
