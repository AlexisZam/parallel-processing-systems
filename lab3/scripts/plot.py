#!/usr/bin/env python3

import pandas as pd
import seaborn as sns

sns.set()

# z1

mapper = {"1": "1st", "2": "2nd"}

data = []
for i in range(3):
    for file in ["z1", "z1-padding"]:
        padding = "padding" in file
        with open(f"../outputs/{file}-{i}.out") as f:
            for line in f:
                if line.startswith("Run"):
                    run = line.split(sep=":")[1].strip()
                elif line.startswith("Nthreads"):
                    tokens = line.split()
                    data.append(
                        [padding, mapper[run], int(tokens[1]), float(tokens[5])]
                    )

columns = ["Padding", "Run", "Number of Threads", "Throughput"]
data = pd.DataFrame(data=data, columns=columns)
sns.relplot(
    x="Number of Threads",
    y="Throughput",
    hue="Run",
    data=data,
    col="Padding",
    kind="line",
    facet_kws={"sharey": False, "legend_out": False},
    marker="o",
).set(xscale="log").set(xticks=[1, 2, 4, 8, 16, 32, 64]).set_xticklabels(
    labels=[1, 2, 4, 8, 16, 32, 64]
).savefig(
    "../plots/z1.pdf", bbox_inches="tight"
)

# z2

mapper = {
    "pthread_lock": "Pthread",
    "tas_lock": "TAS",
    "ttas_lock": "TTAS",
    "array_lock": "Array",
    "clh_lock": "CLH",
}

data = []
for i in range(3):
    with open(f"../outputs/z2-{i}.out") as f:
        for line in f:
            if line.startswith("Lock"):
                lock = line.split(sep=":")[1].strip()
            if line.startswith("List Size"):
                list_size = line.split(sep=":")[1].strip()
            elif line.startswith("Nthreads"):
                if lock != "nosync_lock":
                    tokens = line.split()
                    data.append(
                        [mapper[lock], int(list_size), int(tokens[1]), float(tokens[5])]
                    )

columns = ["Lock", "List Size", "Number of Threads", "Throughput"]
data = pd.DataFrame(data=data, columns=columns)
sns.relplot(
    x="Number of Threads",
    y="Throughput",
    hue="Lock",
    data=data,
    col="List Size",
    col_wrap=2,
    kind="line",
    facet_kws={"sharey": False, "legend_out": False},
    marker="o",
).set(xscale="log").set(xticks=[1, 2, 4, 8, 16, 32, 64]).set_xticklabels(
    labels=[1, 2, 4, 8, 16, 32, 64]
).savefig(
    "../plots/z2.pdf", bbox_inches="tight"
)

# z3

mapper = {
    "ll_fgl": "Fine-Grained",
    "ll_opt": "Optimistic",
    "ll_lazy": "Lazy",
    "ll_nb": "Non-Blocking",
}

data = []
for i in range(2):
    with open(f"../outputs/z3-{i}.out") as f:
        for line in f:
            if line.startswith("LL"):
                ll = line.split(sep=":")[1].strip()
            if line.startswith("List Size"):
                list_size = line.split(sep=":")[1].strip()
            elif line.startswith("Nthreads"):
                tokens = line.split()
                data.append(
                    [
                        mapper[ll],
                        int(list_size),
                        int(tokens[1]),
                        tokens[5],
                        float(tokens[7]),
                    ]
                )

columns = [
    "Synchronization",
    "List Size",
    "Number of Threads",
    "Workload",
    "Throughput",
]
data = pd.DataFrame(data=data, columns=columns)
sns.relplot(
    x="Number of Threads",
    y="Throughput",
    hue="Synchronization",
    data=data,
    row="List Size",
    col="Workload",
    kind="line",
    facet_kws={"sharey": False, "legend_out": False},
    marker="o",
).set(xscale="log").set(xticks=[1, 2, 4, 8, 16, 32, 64]).set_xticklabels(
    labels=[1, 2, 4, 8, 16, 32, 64]
).savefig(
    "../plots/z3.pdf", bbox_inches="tight"
)
