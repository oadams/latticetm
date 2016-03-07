#!/usr/bin/python2
from __future__ import print_function
"""
This is a script for the purpose of more closely examining the transcriptions
and alignments found with latticelm.
"""
import sys
import io
import distance

with io.open(sys.argv[1], encoding="utf8") as plain_f:
    plain_lines = plain_f.readlines()
with io.open(sys.argv[2], encoding="utf8") as tm_f:
    tm_lines = tm_f.readlines()
with io.open(sys.argv[3], encoding="utf8") as gold_f:
    gold_lines = gold_f.readlines()
#with io.open(sys.argv[3], encoding="utf8") as one_best_f:
#    one_best_lines = one_best_f.readlines()

plain_lines = [line.strip() for line in plain_lines]
#one_best_lines = [line.strip() for line in one_best_lines]
tm_lines = [line.strip() for line in tm_lines]
gold_lines = [line.strip() for line in gold_lines]

assert len(plain_lines) == len(tm_lines)
assert len(plain_lines) == len(gold_lines)

i = 0
count_better = 0
count_worse = 0
count_same = 0
for i in range(len(plain_lines)):
    gold_line = gold_lines[i].split()
    plain_line = plain_lines[i].split()
    tm_line = tm_lines[i].split()
    plain_dist = distance.min_edit_distance(gold_line, plain_line)
    tm_dist = distance.min_edit_distance(gold_line, tm_line)

    if tm_dist < plain_dist:
        count_better += 1
    elif tm_dist > plain_dist:
        print(i)
        print("Gold:\t", gold_line)
        print("Plain (med %d):\t" % plain_dist, plain_line)
        #print("1best:\t", one_best_lines[i])
        print("TM (med %d):\t" % tm_dist, tm_line)
        print("===")
        raw_input()
        count_worse += 1
    else:
        count_same += 1
    i += 1

print(count_better)
print(count_worse)
print(count_same)
