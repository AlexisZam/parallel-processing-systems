#!/usr/bin/env python3

import pandas as pd
import seaborn as sns

sns.set()

# z1

data = []
for file in ["naive", "coalesced_A", "reduced_global"]:
    with open(f"../outputs/z1/{file}.out") as f:
        for line in f:
            if line.startswith("/usr/local/cuda/bin/nvcc -D__FLOAT_VALUES"):
                thread_block_x = thread_block_y = 32
                if "-DTHREAD_BLOCK_X=" in line:
                    thread_block_x = int(line.split("-DTHREAD_BLOCK_X=")[1].split()[0])
                    thread_block_y = int(line.split("-DTHREAD_BLOCK_Y=")[1].split()[0])
            elif line.startswith("GPU kernel version"):
                gpu_kernel_version = line.split(":")[1].strip()
            elif line.startswith("Performance"):
                performance = float(line.split(":")[1].strip().split()[0])
                data.append(
                    (
                        gpu_kernel_version,
                        thread_block_x,
                        thread_block_x * thread_block_y,
                        performance,
                    )
                )

columns = ["Kernel", "Thread Block X", "Threads per Block", "Performance"]
data = pd.DataFrame(data=data, columns=columns)

fg = sns.relplot(
    x="Thread Block X",
    y="Performance",
    hue="Threads per Block",
    data=data,
    col="Kernel",
    col_wrap=2,
    legend="full",
    kind="line",
    facet_kws={"sharex": False, "sharey": False, "legend_out": False},
    marker="o",
)
for ax in fg.axes:
    ax.set_xscale("log", basex=2)
fg.savefig("../plots/z1.pdf", bbox_inches="tight")

# z2

data = []
for file in ["naive_1024_1", "coalesced_A_8_16", "reduced_global_32_32", "cublas"]:
    with open(f"../outputs/z2/{file}.out") as f:
        for line in f:
            if line.startswith("Dimension M"):
                M = int(line.split(":")[1].strip())
            elif line.startswith("Dimension N"):
                N = int(line.split(":")[1].strip())
            elif line.startswith("Dimension K"):
                K = int(line.split(":")[1].strip())
            elif line.startswith("GPU kernel version"):
                gpu_kernel_version = line.split(":")[1].strip()
            elif line.startswith("Performance"):
                performance = float(line.split(":")[1].strip().split()[0])
                data.append((gpu_kernel_version, M * N * K, performance))

columns = ["Kernel", "MxNxK", "Performance"]
data = pd.DataFrame(data=data, columns=columns)

fg = sns.relplot(
    x="MxNxK",
    y="Performance",
    hue="Kernel",
    data=data[data["Kernel"] != "cublas"],
    kind="line",
    facet_kws={"legend_out": False},
    marker="o",
)
fg.ax.set_xscale("log", basex=2)
fg.savefig("../plots/z2.pdf", bbox_inches="tight")

fg = sns.relplot(
    x="MxNxK",
    y="Performance",
    hue="Kernel",
    data=data,
    kind="line",
    facet_kws={"legend_out": False},
    marker="o",
)
fg.ax.set_xscale("log", basex=2)
fg.savefig("../plots/z2_cublas.pdf", bbox_inches="tight")
