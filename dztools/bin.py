"""The functions to be used to run infretis via the terminal."""
import pathlib
import sys

from dztools.misc import dzlog, log, get_mapper

MOD_PATH = str(pathlib.Path(__file__).parent.resolve())
FOLDERS = ["funcs", "md"]
MAPPER = get_mapper(FOLDERS, MOD_PATH)
MAPPER["log"] = log
sys.dont_write_bytecode = True

# NOTE: when defining new functionality
# put the import statements in the function defenition
# as to avoid importing loads of libraries, which slows
# down the `dz` call from the command line


def dztool():
    """Map an dztool command to a python function.

    Usage from the command line:
        dz `tool_name` `arguments`

    `tool_name`: the name of the tool you  want to use
    `arguments`: the arguments passed to the tool

    To get more information about each tool, pass `-h` after
    specifying the tool name, e.g. 'inft wham -h'.

    Example usage:
        dz plot -h

    """
    if len(sys.argv) == 1 or sys.argv[1] in ["-h", "help", "--help"]:
        print(dztool.__doc__)
        print("Available inft commands:")
        for key in MAPPER.keys():
            print(f"\t{key}")
        return

    command = "dz " + " ".join(sys.argv[1:]) + "\n"
    dzlog(command, MOD_PATH)
    tool_name = sys.argv[1]
    arguments = sys.argv[2:]

    if tool_name == 'log':
        log(MOD_PATH)
        return

    if tool_name not in list(MAPPER.keys()):
        msg = f"No tool named '{tool_name}', maybe you spelled it wrong?\n \
                \nFor usage information run\n\tinft -h"
        return msg

    tool = MAPPER[tool_name]
    # run the tool function
    tool(arguments)
