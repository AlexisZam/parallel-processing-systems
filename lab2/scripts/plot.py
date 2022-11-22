#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# parse

data = []
for i in range(1, 4):
    with open(f"../outputs/run-{i}.out") as f:
        indices = [0, 2, 4, 6, 8, 12, 14]
        data += [[line.split()[i] for i in indices] for line in f]

columns = ["Method", "X", "Y", "Px", "Py", "Computation Time", "Total Time"]
dtype = dict(zip(columns, [str, str, str, int, int, float, float]))
df = pd.DataFrame(data=data, columns=columns).astype(dtype)

df.replace(
    to_replace={"GaussSeidelSOR": "Gauss-Seidel SOR", "RedBlackSOR": "Red-Black SOR"},
    inplace=True,
)
df["Shape"] = df["X"] + "x" + df["Y"]
df["Number of Processes"] = df["Px"] * df["Py"]
df.drop(columns=["X", "Y", "Px", "Py"], inplace=True)

df = df.groupby(
    by=["Method", "Shape", "Number of Processes"], as_index=False, sort=False
).mean()

# speedup

speedup = []
for _, grouped in df.groupby(by=["Method", "Shape"], as_index=False, sort=False):
    speedup += list(
        grouped[grouped["Number of Processes"] == 1]["Total Time"].iloc[0]
        / grouped["Total Time"]
    )
df["Speedup"] = speedup

sns.set(style="whitegrid")

for shape, data in df.groupby(by="Shape"):
    sns.lineplot(
        x="Number of Processes", y="Speedup", hue="Method", data=data, marker="o"
    )
    plt.title(f"Shape = {shape}")
    plt.savefig(f"../plots/speedup-{shape}.pdf", bbox_inches="tight")
    plt.close()

# time

df = df[df["Number of Processes"] >= 8]
df.replace(
    to_replace={
        "Gauss-Seidel SOR": "Gauss-Seidel\nSOR",
        "Red-Black SOR": "Red-Black\nSOR",
    },
    inplace=True,
)

df = df.melt(
    id_vars=["Method", "Shape", "Number of Processes"],
    value_vars=["Total Time", "Computation Time"],
    var_name="Time Type",
    value_name="Time",
)

for shape, data in df.groupby(by="Shape"):
    sns.catplot(
        x="Method",
        y="Time",
        hue="Time Type",
        data=data,
        col="Number of Processes",
        kind="bar",
        aspect=0.7,
    )
    plt.savefig(f"../plots/time-{shape}.pdf", bbox_inches="tight")
    plt.close()
