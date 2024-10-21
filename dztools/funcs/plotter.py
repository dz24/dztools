from typing import Annotated

import typer

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


def plot(
    i: Annotated[list[str], typer.Option("-i", help="input")],
    log: Annotated[bool, typer.Option("-log", help="input")] = False,
    sci: Annotated[bool, typer.Option("-sci", help="input")] = False,
    avg: Annotated[bool, typer.Option("-avg", help="input")] = False,
    hl: Annotated[
        list[float], typer.Option("-hl", help="horizontal lines")
    ] = [],
    dim: Annotated[list[int], typer.Option("-dim", help="number of.")] = [1],
):
    """
    Plot via the terminal using matplotlib.
    """
    import matplotlib.pyplot as plt
    import matplotlib
    # matplotlib.use('module://matplotlib-backend-wezterm')
    import numpy as np

    for idx, file in enumerate(i):
        data = np.loadtxt(file)
        alpha = 1 if not avg else 0.2
        for dim0 in dim:
            plt.plot(data[:, 0], data[:, dim0], alpha=alpha)
        if avg:
            plt.axhline(np.average(data[:, 1]), color=COLS[idx])
            print(i, np.average(data[:, 1]))
    if log:
        plt.yscale("log")
    if hl is not None:
        for hline in hl:
            plt.axhline(hline, color="k")

    if sci:
        import scienceplots
        plt.style.use("science")
        plt.savefig("fig.pdf")
    else:
        plt.show()
