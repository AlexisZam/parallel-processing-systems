#!/usr/bin/env python3

import pandas as pd

data = []
for i in range(1, 4):
    with open(f"../outputs/run-conv-{i}.out") as f:
        indices = [0, 12, 14, 16]
        data += [[line.split()[i] for i in indices] for line in f]

columns = ["Method", "Computation Time", "Total Time", "Convergence Time"]
dtype = dict(zip(columns, [str, float, float, float]))
df = pd.DataFrame(data=data, columns=columns).astype(dtype)

df.replace(
    to_replace={"GaussSeidelSOR": "Gauss-Seidel SOR", "RedBlackSOR": "Red-Black SOR"},
    inplace=True,
)
df = df[["Method", "Total Time", "Computation Time", "Convergence Time"]]

df = df.groupby(by="Method", sort=False).mean().round(decimals=6)

df.to_csv(path_or_buf="../tables/run-conv.csv")
