import glob
import importlib
from inspect import getmembers, isfunction
from os.path import basename, isfile

MAXLINES = 100


def get_mapper(folders, mod_path):
    mapper = {}
    for folder in folders:
        fpath = mod_path + f"/{folder}"
        files = [basename(i)[:-3] for i in glob.glob(fpath + "/*py")]
        for file in files:
            mod = importlib.import_module(
                f".{folder}.{file}", package="dztools"
            )
            funcs = getmembers(mod, isfunction)
            for func in funcs:
                mapper[func[0]] = func[1]
    return mapper


def dzlog(command, mod_path):
    commands = [command]
    log_path = mod_path + "/.log"
    if isfile(log_path):
        with open(log_path, "r") as read:
            for line in read:
                commands.append(line)
    commands = list(set(commands))[-MAXLINES:]

    with open(log_path, "w") as write:
        for line in commands:
            write.write(line)


def log(mod_path):
    log_path = mod_path + "/.log"
    with open(log_path, 'r') as read:
        for line in read:
            print(line.rstrip())
