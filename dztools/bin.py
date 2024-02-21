"""The functions to be used to run infretis via the terminal."""
import sys

from dztools.calc_dens import (
    calc_dens_box,
)
from dztools.center_xyz import (
    center_periodic,
)
from dztools.log import (
    dzlog,
    dzlog_print,
)
from dztools.orderp import (
    count_w,
    co2_op,
    calc_dist,
)
from dztools.plotter import(
    plot_data,
)


# NOTE: when defining new functionality
# put the import statements in the function defenition
# as to avoid importing loads of libraries, which slows
# down the `inft` call from the command line

MAPPER = {
        "dens": calc_dens_box,
        "pbc": center_periodic,
        "log": dzlog_print,
        "cw": count_w,
        "co2_op": co2_op,
        "calc_dist": calc_dist,
        "plot": plot_data,
}

def dztool():
    """Map an dztool command to a python function.

    Usage from the command line:
        inft `tool_name` `arguments`

    `tool_name`: the name of the tool you  want to use
    `arguments`: the arguments passed to the tool

    To get more information about each tool, pass `-h` after
    specifying the tool name, e.g. 'inft wham -h'.

    Example usage:
        inft wham -toml infretis.toml -data infretis_data.txt
        inft concatenate -h

    """
    if len(sys.argv) == 1 or sys.argv[1] in ["-h","help","--help"]:
        print(dztool.__doc__)
        print("Available inft commands:")
        for key in MAPPER.keys():
            print(f"\t{key}")
        return

    command = 'dz ' + ' '.join(sys.argv[1:]) + '\n'
    dzlog(command)
    tool_name = sys.argv[1]
    arguments = sys.argv[2:]

    if tool_name not in list(MAPPER.keys()):
        msg = f"No tool named '{tool_name}', maybe you spelled it wrong?\n \
                \nFor usage information run\n\tinft -h"
        return msg

    tool = MAPPER[tool_name]
    # run the tool function
    tool(arguments)
