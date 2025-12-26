"""The functions to be used to run dztools via the terminal."""

import pathlib
import sys

import typer

from dztools.misc.bin_help import dzlog, get_mapper, log

# avoids the automatic func_name -> func-name conversion.
# https://github.com/fastapi/typer/issues/341#issuecomment-1687580973
typer.main.get_command_name = lambda name: name

# define constants
MOD_PATH = str(pathlib.Path(__file__).parent.resolve())
# FOLDERS = ["funcs", "md", "mem", "chem", "ruru"]
FOLDERS = ["funcs", "md", "mem", "chem", "ruru", "co2"]
MAPPER = get_mapper(FOLDERS, MOD_PATH)

# create the typer app
app = typer.Typer(
    no_args_is_help=True,
    help="dz24's dztools CLI",
    epilog="katt",
    context_settings={"help_option_names": ["-h", "--help"]},
)

# decorating imported mapper functions
for func in MAPPER.values():
    app.command()(func)

# NOTE: when defining new functionality
# put the import statements in the function defenition
# as to avoid importing loads of libraries, which slows
# down the `dz` call from the command line
