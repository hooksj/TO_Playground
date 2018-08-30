from sys import platform
import subprocess
import os
import argparse

# for the timeout thing at the end
if platform == "linux" or platform == "linux2":
    # linux
    import sys
    from select import select


    def wait_for_input(timeout):
        print "Press ENTER to continue or wait " + str(timeout) + " seconds..."
        rlist, _, _ = select([sys.stdin], [], [], timeout)

        if rlist:
            print "Continuing..."
        else:
            print "Timed out..."

elif platform == "darwin":
    # OS X
    def wait_for_input(timeout):
        print "function isn't configured for osx"

elif platform == "win32":
    # Windows...
    import sys, time, msvcrt


    def wait_for_input(timeout):
        startTime = time.time()
        inp = None
        print "Press ENTER to continue or wait " + str(timeout) + " seconds..."
        while True:
            if msvcrt.kbhit():
                inp = msvcrt.getch()
                break
            elif time.time() - startTime > timeout:
                break

        if inp:
            print "Config selected..."
        else:
            print "Timed out..."


def check_limit(val, limit):
    """
    Limits values
    :param val:value to limit
    :param limit: tuple (lower_limit, upper_limit)
    :return: the limited value
    """
    # print "check limit"
    if val < limit[0]:
        return limit[0]
    if val > limit[1]:
        return limit[1]
    else:
        return val


# def play_sound(filepath):
#     Process(target=playsound, args=(filepath))

def play_mp3(path):
    subprocess.Popen(['mpg123', '-q', path])


def play_glass_ping():
    util_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = util_dir + "/../data/mp3/glass_ping-Go445-1207030150.mp3"
    play_mp3(filepath)


def play_beep():
    util_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = util_dir + "/../data/mp3/beep.mp3"
    play_mp3(filepath)


def player_arg_parser(filename):
    parser = argparse.ArgumentParser(description=filename)
    parser.add_argument('-s', '--simulation',
                        help='run a VREP simulation',
                        action='store_true')
    parser.add_argument('-d', '--dynamixel',
                        help='use Dynamixels',
                        action='store_true')
    args = parser.parse_args()
    if not args.dynamixel and not args.simulation:
        print(
            """
            ERROR: You need to specify either a simulation or a robot to output your keyframes to."
        
            usage: python {} [-h] [-s] [-d]  
        
            optional arguments:
            -h, --help        show this help message and exit
            -s, --simulation  run a VREP simulation
            -d, --dynamixel   use Dynamixels
            """.format(filename))
        quit()

    if args.dynamixel and args.simulation:
        print("ERROR: Only select either dynamixel or simulation for now.")
        quit()

    return args