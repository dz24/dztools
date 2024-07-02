import typer
from typing_extensions import Annotated
# app = typer.Typer()

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

# @app.command()
def plot(
    i: Annotated[list[str], typer.Option(help="input")],
    dim: Annotated[int, typer.Option(help="number of.")] = 1,
    log: bool = False,
    sci: bool = False,
    avg: bool = False,
    hl: list[float] = []
):
    """
    Plot via the terminal using matplotlib.
    """
# def plot(arguments):
    # parser = argparse.ArgumentParser(
    #     description="Make xyz file periodic for Jmol viewing."
    # )

    # parser.add_argument(
    #     "-i",
    #     help="Input file",
    #     nargs="+",
    # )

    # parser.add_argument(
    #     "-dim", help="y dims", type=int, nargs="+", default=[1]
    # )

    # parser.add_argument(
    #     "-log",
    #     help="Log scale",
    #     action="store_true",
    # )

    # parser.add_argument(
    #     "-sci",
    #     help="science plot",
    #     action="store_true",
    # )

    # parser.add_argument(
    #     "-avg",
    #     help="avg",
    #     action="store_true",
    # )
    # parser.add_argument(
    #     "-hl",
    #     help="horizontal lines",
    #     type=float,
    #     nargs="+",
    #     required=False,
    # )

    # args = parser.parse_args(arguments)
    # print('burm')
    # data = np.loadtxt(i)
    # return

    if sci:
        import scienceplots
        plt.style.use("science")

    # data = np.loadtxt(i)
    # plt.plot(data[:, 0], data[:, 1])
    for file in i:
        data = np.loadtxt(file)
        plt.plot(data[:, 0], data[:, dim])
        if avg:
            print(i, np.average(data[:, 1]))
    # for dim0 in dim:
    # k    plt.plot(data[:, 0], data[:, dim0])
    if log:
        plt.yscale("log")
    if hl is not None:
        for hline in hl:
            plt.axhline(hline, color="k")

    if sci:
        plt.savefig("fig.pdf", bbox_inches="tight")
    else:
        plt.show()
