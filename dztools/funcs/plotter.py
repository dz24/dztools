import argparse

import matplotlib.pyplot as plt
import numpy as np

COLS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]


def plot(arguments):
    parser = argparse.ArgumentParser(
        description="Make xyz file periodic for Jmol viewing."
    )

    parser.add_argument(
        "-i",
        help="Input file",
        nargs="+",
    )

    parser.add_argument(
        "-dim", help="y dims", type=int, nargs="+", default=[1]
    )

    parser.add_argument(
        "-log",
        help="Log scale",
        action="store_true",
    )

    parser.add_argument(
        "-sci",
        help="science plot",
        action="store_true",
    )

    parser.add_argument(
        "-avg",
        help="avg",
        action="store_true",
    )
    parser.add_argument(
        "-hl",
        help="horizontal lines",
        type=float,
        nargs="+",
        required=False,
    )

    args = parser.parse_args(arguments)

    if args.sci:
        import scienceplots
        plt.style.use("science")

    for i in args.i:
        data = np.loadtxt(i)
        if args.avg:
            print(i, np.average(data[:, 1]))
        for dim in args.dim:
            plt.plot(data[:, 0], data[:, dim])
    if args.log:
        plt.yscale("log")
    if args.hl is not None:
        for hline in args.hl:
            plt.axhline(hline, color="k")

    if args.sci:
        plt.savefig("fig.pdf", bbox_inches="tight")
    else:
        plt.show()