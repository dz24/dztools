"""The functions to be used to run dztools via the terminal."""

import pathlib
import sys

import typer

from dztools.misc.bin_help import dzlog, get_mapper, log

# define constants
MOD_PATH = str(pathlib.Path(__file__).parent.resolve())
FOLDERS = ["funcs", "md"]
MAPPER = get_mapper(FOLDERS, MOD_PATH)
MAPPER["log"] = log

# log the command
dzlog(" ".join(sys.argv) + "\n", MOD_PATH)

# create the typer app
app = typer.Typer(
    no_args_is_help=True,
    help="ape",
    epilog="katt",
    context_settings={"help_option_names": ["-h", "--help"]},
)

# decorating imported mapper functions
for func in MAPPER.values():
    app.command()(func)

# manually add function from other folder
app.command()(log)

# NOTE: when defining new functionality
# put the import statements in the function defenition
# as to avoid importing loads of libraries, which slows
# down the `dz` call from the command line
