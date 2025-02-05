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
    i: Annotated[list[str], typer.Option("-i")],
    log: Annotated[bool, typer.Option("-log")] = False,
    sci: Annotated[bool, typer.Option("-sci")] = False,
    avg: Annotated[bool, typer.Option("-avg")] = False,
    leg: Annotated[bool, typer.Option("-leg")] = False,
    hl: Annotated[
        list[float], typer.Option("-hl", help="horizontal lines")
    ] = [],
    vl: Annotated[
        list[float], typer.Option("-vl", help="vertical lines")
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

    print("ww")

    if sci:
        import scienceplots
        plt.style.use("science")

    for idx, file in enumerate(i):
        data = np.loadtxt(file, comments=["@", "#"])
        alpha = 1 if not avg else 0.2
        for dim0 in dim:
            plt.plot(data[:, 0], data[:, dim0], alpha=alpha, label=str(idx))
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

def plot_ps(
    pns: Annotated[list[str], typer.Option("-pns", help="input")],
    log: Annotated[bool, typer.Option("-log", help="input")] = False,
    sci: Annotated[bool, typer.Option("-sci", help="input")] = False,
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
    import os
    import numpy as np

    if sci:
        import scienceplots
        plt.style.use("science")

    for idx, pn in enumerate(pns):
        print(pn)
        print('tiger', os.path.join(pn, "order.txt"))
        pn_o = np.loadtxt(os.path.join(pn, "order.txt"))
        pn_t = np.loadtxt(os.path.join(pn, "traj.txt"), dtype=str)
        trajs = list(set(pn_t[:, 1]))
        for idx0, traj in enumerate(trajs):
            alpha = 0.5 if traj in (pn_t[0, 1], pn_t[-1, 1]) else 1
            print(alpha, idx0, len(trajs)-1)
            x= pn_o[:, 0][pn_t[:, 1] == traj]
            y= pn_o[:, 1][pn_t[:, 1] == traj]
            if len(x) > 1:
                plt.plot(x*0.01, y, color=f"C{idx0}", alpha=alpha)
            else:
                plt.plot(np.array([x[0]-1, x[0], x[0]+1])*0.01, [pn_o[:, 1][int(x[0]-1)],y[0],pn_o[:, 1][int(x[0]+1)]], color=f"C{idx0}", alpha=alpha)
        plt.xlabel("Time [ns]")
        plt.ylabel("Chain Order Parameter")
        if sci:
            plt.savefig("fig.pdf")
        else:
            plt.show()

        print(trajs)
        exit()
        for dim0 in dim:
            plt.plot(data[:, 0], data[:, dim0], alpha=alpha)
        if avg:
            plt.axhline(np.average(data[:, 1]), color=COLS[0])
            print(i, np.average(data[:, 1]))
    if log:
        plt.yscale("log")
    if hl is not None:
        for hline in hl:
            plt.axhline(hline, color="k")

    if sci:
        plt.savefig("fig.pdf")
    else:
        plt.show()
