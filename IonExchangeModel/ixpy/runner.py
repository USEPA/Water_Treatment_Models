# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 08:20:32 2020

Run script for ixpy models

@author: LHaupert
"""

import argparse
from .hsdmix import HSDMIX
from .psdmix import PSDMIX


def parse_args(args):
    """
    Parse arguments from command line.

    Returns
    -------
    Argument object to pass to run_HSDMIX
    """
    parser = argparse.ArgumentParser(prog='ixpy')
    parser.add_argument('mode', choices=['hsdm', 'psdm'])
    parser.add_argument('input_fname')
    parser.add_argument('output_fname')
    parser.add_argument('-t', '--t_unit')
    parser.add_argument('-c', '--c_unit')
    parsed_args = parser.parse_args(args)
    return parsed_args


def run_ixpy(args):
    input_fname = args.input_fname
    output_fname = args.output_fname
    option_dict = {}
    if args.t_unit:
        option_dict.update({'period':args.t_unit})
    if args.c_unit:
        option_dict.update({'units':args.c_unit})
        
    if args.mode == 'hsdm':
        IEX = HSDMIX(input_fname)
    elif args.mode == 'psdm':
        IEX = PSDMIX(input_fname)

    IEX.solve()
    IEX.save_results(output_fname, **option_dict)

    return None
    
