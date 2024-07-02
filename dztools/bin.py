"""The functions to be used to run infretis via the terminal."""
import pathlib
import sys
import typer
# app = typer.Typer()

from dztools.misc.bin_help import dzlog, log, get_mapper

MOD_PATH = str(pathlib.Path(__file__).parent.resolve())
FOLDERS = ["funcs", "md"]
MAPPER = get_mapper(FOLDERS, MOD_PATH)
# app = typer.Typer()
# import dztools.funcs.plotter
# app.add_typer(dztools.funcs.plotter.app, name="plot")

# @app.command()
# def shoot():
#     """
#     Shoot the portal gun
#     """
#     typer.echo("Shooting portal gun")
# @app.command()
# def plot():
#     """
#     """
#     typer.run(MAPPER['plot'])
    # typer.echo("Shooting portal gun")
# for name, func in MAPPER.items():
#     app.add_typer(func, name=name)

# print(MAPPER)
# exit('f')
MAPPER["log"] = log

# NOTE: when defining new functionality
# put the import statements in the function defenition
# as to avoid importing loads of libraries, which slows
# down the `dz` call from the command line

# @app.command()
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
    # app()
    # return
    if len(sys.argv) == 1 or sys.argv[1] in ["-h", "help", "--help"]:
        print(dztool.__doc__)
        print("Available inft commands:")
        for key in MAPPER.keys():
            print(f"\t{key}")
        return

    command = "dz " + " ".join(sys.argv[1:]) + "\n"
    dzlog(command, MOD_PATH)
    tool_name = sys.argv[1]
    # arguments = sys.argv[2:]
    # print('pear', sys.argv)
    sys.argv.pop(1)

    if tool_name == 'log':
        log(MOD_PATH)
        return

    if tool_name not in list(MAPPER.keys()):
        msg = f"No tool named '{tool_name}', maybe you spelled it wrong?\n \
                \nFor usage information run\n\tinft -h"
        return msg

    tool = MAPPER[tool_name]
    # run the tool function
    typer.run(tool)
    # tool(arguments)
    # tool(arguments)
