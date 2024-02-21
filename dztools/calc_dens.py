import argparse

# dens =  g / V

MMASS_W = 18.01528  # g / mol water
AVO = 6.0221408e+23

def calc_dens_box(arguments):
    parser = argparse.ArgumentParser(
        description="Calculate the theta and phi angle for an \
                .sdf file given a set of indices"
    )

    parser.add_argument(
        "-num_w", help="The number of water molecules"
    )
    parser.add_argument(
        "-c",
        help="The length of a cubic cell",
        type=int,
        nargs="+",
    )

    args = parser.parse_args(arguments)
    print(args)
    print(arguments)

    # mmass_w = 18.01528  # g / mol water
    # n_w = 64

    # mass_w = (n_w/avo)*mmass_w
    # # cell = 13.694596613173969  # angstrom
    # cell = 12.4138
    # ang_to_cm = 10**(-8)
    # vol = (cell*ang_to_cm)**3 #cm3
    # dens = mass_w/vol
    # print('whada')
    # print('mass:', mass_w, 'g')
    # print('vol :', vol, 'cm3')
    # print('dens:', dens, 'g/cm3')
