# -*- coding: utf-8 -*-

from .hsdmix import run_HSDMIX, parse_args

def main():
    args = parse_args()
    run_HSDMIX(args)

if __name__=='__main__':
    main()
    
