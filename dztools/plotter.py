import numpy as np
import matplotlib.pyplot as plt
import argparse


def plot_data(arguments):
    parser = argparse.ArgumentParser(
        description="Make xyz file periodic for Jmol viewing."
    )

    parser.add_argument(
        "-i",
        help="Input file",
        nargs="+",
    )

    parser.add_argument(
        "-log",
        help="Log scale",
        action='store_true',
    )

    args = parser.parse_args(arguments)
    for i in args.i:
        data = np.loadtxt(i)
        plt.plot(data[:, 0], data[:, 1])
        # plt.plot(data[:, 0], data[:, 3])
    # plt.plot(data[:, 0], data[:, 2])
    # plt.plot(data[:, 0], data[:, 3])
    if args.log:
        plt.yscale('log')
    plt.show()

