import pathlib

MOD_PATH = str(pathlib.Path(__file__).parent.resolve())
LOG_PATH = MOD_PATH + "/funcs/.log"
MAXLINES = 100


def dzlog(command):
    commands = [command]
    with open(LOG_PATH, 'r') as read:
        for line in read:
            commands.append(line)
    commands = set(commands)

    with open(LOG_PATH, 'w') as write:
        for line in commands:
            write.write(line)


# def dzlog_print(command):
#     parser = argparse.ArgumentParser(
#         description="Make xyz file periodic for Jmol viewing."
#     )
#     with open(LOG_PATH, 'r') as read:
#         for line in read:
#             print(line.rstrip())


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
