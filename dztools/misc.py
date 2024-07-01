from os.path import basename as bs
import glob
import importlib
from inspect import getmembers, isfunction

MAXLINES = 100

def get_mapper(folders, mod_path):
    mapper = {}
    for folder in folders:
        fpath = mod_path + f"/{folder}"
        files = [bs(i)[:-3] for i in glob.glob(fpath + "/*py")]
        for file in files:
            mod = importlib.import_module(f'.{folder}.{file}',
                                          package='dztools')
            funcs = getmembers(mod, isfunction)
            for func in funcs:
                mapper[func[0]] =  func[1]
    return mapper


def dzlog(command, mod_path):
    commands = [command]
    log_path = mod_path + "/.log"
    with open(log_path, 'a+') as read:
        for line in read:
            commands.append(line)
    commands = list(set(commands))[-MAXLINES:]

    with open(log_path, 'w') as write:
        for line in commands:
            write.write(line)


#         "dens": calc_dens_box,
#         "pbc": center_periodic,
#         "log": dzlog_print,
#         "cw": count_w,
#         "co2_op": co2_op,
#         "calc_dist": calc_dist,
#         "plot": plot_data,
#         "ruru_op": calc_op_ruru,
#         "ruru_dist": calc_ruru_dist,
#         "ruru_cat": ruru_cat,
#         "ruru_x": calc_ruru_x,
#         "ruru_idx": find_rel_idxes,

