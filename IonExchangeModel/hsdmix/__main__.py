# -*- coding: utf-8 -*-

import sys

from .runner import run_HSDMIX, parse_args

def main():
    args = parse_args(sys.argv[1:])
    run_HSDMIX(args)

if __name__=='__main__':
    main()
    
