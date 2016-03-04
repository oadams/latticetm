#!/usr/bin/python3
"""
This is a script for the purpose of more closely examining the transcriptions
and alignments found with latticelm.
"""
import sys
import io

with io.open(sys.argv[1], encoding="utf8") as plain_f:
    plain_lines = plain_f.readlines()
with io.open(sys.argv[2], encoding="utf8") as tm_one_best_f:
    tm_one_best_lines = tm_one_best_f.readlines()
#with io.open(sys.argv[3], encoding="utf8") as one_best_f:
#    one_best_lines = one_best_f.readlines()

plain_lines = [line.strip() for line in plain_lines]
#one_best_lines = [line.strip() for line in one_best_lines]
tm_one_best_lines = [line.strip() for line in tm_one_best_lines]

assert len(plain_lines) == len(tm_one_best_lines)
#assert len(plain_lines) == len(one_best_lines)
i = 0
count = 0 
for i in range(len(plain_lines)):
    if plain_lines[i].strip() != tm_one_best_lines[i]:
        print("Sentence %d is different" % i)
        print("Plain:\t", plain_lines[i])
        #print("1best:\t", one_best_lines[i])
        print("TM:\t", tm_one_best_lines[i])
        print("===")
        count += 1
        input()
print(count, "/", len(plain_lines))
