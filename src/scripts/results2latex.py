#!/usr/bin/python3

import sys

results_path = sys.argv[1]
with open(results_path) as results_f:
    lines = results_f.readlines()

results = []
i = 2
best = 1
best_is = []
while i < len(lines):
    result = []
    result.append(lines[i])
    i += 1
    result.extend(lines[i].strip().split())
    result[0] = float("%.3f" % float(result[0]))
    if result[0] < best:
        best = result[0]
        best_is = [len(results)]
    if result[0] == best:
        best_is.append(len(results))
    i += 1
    results.append(result)

asr = float(lines[0].split()[2])
print("\\begin{tabular}{| c c | c |}")
print("\\hline")
print("$\\alpha$ & $\\lambda$ & WER\\\\")
print("\\hline")
i = 0
for result in results:
    if i in best_is:
        print("%s & %s & \\textbf{%.3f}\\\\" % (result[2], result[3], result[0]))
    else:
        print("%s & %s & %.3f\\\\" % (result[2], result[3], result[0]))
    i += 1
print("\\hline")
print("\\multicolumn{2}{| r |}{ASR 1-best} & %.3f\\\\" % asr)
print("\\hline")
print("\\end{tabular}")
