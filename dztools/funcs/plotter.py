from typing import Annotated as And
from typer import Option as Opt

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
    i: And[list[str], Opt("-i")],
    log: And[bool, Opt("-log")] = False,
    sci: And[bool, Opt("-sci")] = False,
    avg: And[bool, Opt("-avg")] = False,
    leg: And[bool, Opt("-leg")] = False,
    scat: And[bool, Opt("-scat")] = False,
    hl: And[ list[float], Opt("-hl", help="horizontal lines") ] = [],
    vl: And[ list[float], Opt("-vl", help="vertical lines") ] = [],
    dim: And[list[int], Opt("-dim", help="number of.")] = [1],
):
    """
    Plot via the terminal using matplotlib.
    """
    import matplotlib.pyplot as plt
    import matplotlib
    # matplotlib.use('module://matplotlib-backend-wezterm')
    import numpy as np

    if sci:
        import scienceplots
        plt.style.use("science")

    for idx, file in enumerate(i):
        data = np.loadtxt(file, comments=["@", "#"])
        alpha = 1 if not avg else 0.2
        for dim0 in dim:
            plt.plot(data[:, 0], data[:, dim0], alpha=alpha, label=str(idx))
            if scat:
                plt.scatter(data[:, 0], data[:, dim0], alpha=alpha, label=str(idx))
        if avg:
            plt.axhline(np.average(data[:, 1]), color=COLS[idx])
            print(i, np.average(data[:, 1]))
    if log:
        plt.yscale("log")
    if hl is not None:
        for hline in hl:
            plt.axhline(hline, color="k")
    if vl is not None:
        for vline in vl:
            plt.axvline(vline, color="k")

    if leg:
        plt.legend(frameon=False)
    if sci:
        plt.savefig("fig.pdf")
    else:
        plt.show()
