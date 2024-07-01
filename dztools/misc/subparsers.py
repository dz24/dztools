import argparse

def add_io(parser, addin=True, addout=True):
    print('w')
    if addin:
        parser.add_argument(
            "-i", help="Input file"
        )
    if addout:
        parser.add_argument(
            "-o",
            help="Output file",
            nargs='?',
            const='periodic.xyz',
        )
