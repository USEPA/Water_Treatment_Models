# -*- coding: utf-8 -*-

import sys

from .runner import parse_args, run_ixpy

def main():
    args = parse_args(sys.argv[1:])
    run_ixpy(args)

if __name__=='__main__':
    main()
    
