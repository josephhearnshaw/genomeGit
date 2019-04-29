#!/usr/bin/env python2
import subprocess
import sys


def check_mod(module):
    """
    Check to see if the module is present within installed modules

    Requires the module

    """
    if module in installed_mod:
        pass
    else:
        print("\n***{} wasn't found for Python 2.7+,"
              " you may need to install it for "
              "genomeGit 3.0 to work properly.***\n".format(module))


# Obtain output form pip as list
required_pkg = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
# split it up into seperate strings
installed_mod = [pkg.decode().split('==')[0] for pkg in required_pkg.split()]
# define modules required
modules = ['pytabix', 'pyfaidx']
# iterate through the modules and  check them.
for package in modules:
    check_mod(package)
