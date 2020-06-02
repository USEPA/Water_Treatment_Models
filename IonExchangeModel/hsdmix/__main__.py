# -*- coding: utf-8 -*-

from .hsdmix import run_HSDMIX, parse_args

if __name__=='__main__':
    args = parse_args()
    run_HSDMIX(args)
    
