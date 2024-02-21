import pathlib

MOD_PATH = str(pathlib.Path(__file__).parent.resolve())
LOG_PATH = MOD_PATH + "/log/log.txt"
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
    

def dzlog_print(arguments):
    with open(LOG_PATH, 'r') as read:
        for line in read:
            print(line.rstrip())

    
