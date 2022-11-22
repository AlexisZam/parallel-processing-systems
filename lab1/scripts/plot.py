#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os.path
import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: ./{} source destination".format(os.path.basename(sys.argv[0])))
    exit(1)

tuples = []
serialTime = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith("Number of threads:"):
            numberOfThreads = int(line.split()[3])
        else:
            # size = int(line.split()[2])
            size = int(line.split(",")[1])
            # time = float(line.split()[6])
            time = float(line.split(",")[-1])
            if numberOfThreads == 1:
                serialTime[size] = time
            speedup = serialTime[size] / time
            tuples.append((size, numberOfThreads, time, speedup))
tuples.sort()

df = pd.DataFrame(tuples, columns=("size", "numberOfThreads", "time", "speedup"))
print(df)

configuration = "-".join(os.path.basename(sys.argv[1]).split(".")[0].split("-")[1:])

for size, group in df.groupby("size"):
    group.plot(
        "numberOfThreads",
        "time",
        legend=False,
        marker="o",
        title="Size: {}".format(size),
    )
    plt.xlabel("Number of threads")
    plt.ylabel("Time")
    if not configuration:
        plt.savefig(os.path.join(sys.argv[2], "time-{}.pdf".format(size)))
    else:
        plt.savefig(
            os.path.join(sys.argv[2], "{}-time-{}.pdf".format(configuration, size))
        )

fig, ax = plt.subplots()
for size, group in df.groupby("size"):
    group.plot("numberOfThreads", "speedup", ax=ax, label=size, marker="o")
plt.legend(title="Size")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
if not configuration:
    plt.savefig(os.path.join(sys.argv[2], "speedup.pdf"))
else:
    plt.savefig(os.path.join(sys.argv[2], "{}-speedup.pdf".format(configuration)))
