#!/usr/bin/python2
from __future__ import print_function
"""
This is a script for the purpose of more closely examining the transcriptions
and alignments found with latticelm.
"""
import sys
import io
import distance

with io.open(sys.argv[1], encoding="utf8") as gold_f:
    gold_lines = gold_f.readlines()
with io.open(sys.argv[2], encoding="utf8") as plain_f:
    plain_lines = plain_f.readlines()
with io.open(sys.argv[3], encoding="utf8") as tm_f:
    tm_lines = tm_f.readlines()
with io.open(sys.argv[4], encoding="utf8") as giza_f:
    giza_lines = giza_f.readlines()

gold_lines = [line.strip() for line in gold_lines]
plain_lines = [line.strip() for line in plain_lines]
tm_lines = [line.strip() for line in tm_lines]
giza_lines = [line.strip() for line in giza_lines]

assert len(plain_lines) == len(tm_lines)
assert len(plain_lines) == len(gold_lines)
assert len(plain_lines) == len(giza_lines)

i = 0
count_better = 0
count_worse = 0
count_same = 0

giza_better = 0
giza_worse = 0
giza_same = 0

both_better = 0
both_worse = 0
giza_better_tm_worse = 0
for i in range(len(plain_lines)):
    gold_line = gold_lines[i].split()
    plain_line = plain_lines[i].split()
    tm_line = tm_lines[i].split()
    giza_line = giza_lines[i].split()

    plain_dist = distance.min_edit_distance(gold_line, plain_line)
    tm_dist = distance.min_edit_distance(gold_line, tm_line)
    giza_dist = distance.min_edit_distance(gold_line, giza_line)

    if giza_dist < plain_dist:
        giza_better += 1
        if tm_dist < plain_dist:
            both_better += 1
        if tm_dist > plain_dist:
            giza_better_tm_worse += 1
            #print(i)
            #print("Gold:\t", gold_line)
            #print("Plain (med %d):\t" % plain_dist, plain_line)
            #print("GIZA (med %d):\t" % giza_dist, giza_line)
            #print("TM (med %d):\t" % tm_dist, tm_line)
            #print("===")
            #raw_input()
    elif giza_dist > plain_dist:
        giza_worse += 1
        if tm_dist > plain_dist:
            both_worse += 1
    else:
        giza_same += 1

    if tm_dist < plain_dist:
        count_better += 1
        if tm_dist < giza_dist:
            print(i)
            print("Gold:\t", gold_line)
            print("Plain (med %d):\t" % plain_dist, plain_line)
            print("GIZA (med %d):\t" % giza_dist, giza_line)
            print("TM (med %d):\t" % tm_dist, tm_line)
            print("===")
            raw_input()
            pass
    elif tm_dist > plain_dist:
        count_worse += 1
    else:
        count_same += 1
    i += 1

print("\tLat-TM\tGIZA\tboth")
print("Better:\t%d\t%d\t%d" % (count_better,giza_better,both_better))
print("Worse:\t%d\t%d\t%d" % (count_worse,giza_worse,both_worse))
print("Same:\t%d\t%d" % (count_same,giza_same))
print("Giza better TM worse: %d" % giza_better_tm_worse)
