#!/usr/bin/env python3

"""
Script used to execute RapidMoc package from command line.

"""


import sys

from rapidmoc.rapidmoc import main

if __name__ == '__main__':
    
    try:
        main()
    except KeyboardInterrupt as err:
        print(err)
        sys.exit()
        
        
