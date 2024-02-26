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

    parser.add_argument(
        "-avg",
        help="avg",
        action='store_true',
    )
    parser.add_argument(
        "-hl",
        help="horizontal lines",
        type=int,
        nargs="+",
        required=False,
    )

    args = parser.parse_args(arguments)
    for i in args.i:
        data = np.loadtxt(i)
        if args.avg:
            print(i, np.average(data[:, 1]))
        plt.plot(data[:, 0], data[:, 1])
        # plt.plot(data[:, 0], data[:, 3])
    # plt.plot(data[:, 0], data[:, 2])
    # plt.plot(data[:, 0], data[:, 3])
    if args.log:
        plt.yscale('log')
    if args.hl is not None:
        for hline in args.hl:
            plt.axhline(hline, color='k')

    plt.show()

