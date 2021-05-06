#!/usr/bin/env python3

import subprocess, argparse
parser = argparse.ArgumentParser(
      description='convenience script to compile executables')
parser.add_argument(
      '-d','--debug',
      dest='debug',
      type=str,
      default='',
      help='Whether to add debugging information (true/false)')
parser.add_argument('-p','--profile',
      dest='profile',
      type=str,
      default='',
      help='Whether to add profiling information (true/false)')

args = parser.parse_args()

print("=====================================================")
print("Compiling radial Chebyshev executable")
print("=====================================================")
subprocess.call(
      'make debug={} profile={} use_cheb=true use_contiguous_longitudes=true'.format(args.debug, args.profile), 
      shell=True
      )
print("=====================================================")
print("Compiling radial finite difference executable")
print("=====================================================")
subprocess.call(
      'make debug={} profile={} use_contiguous_longitudes=true'.format(args.debug, args.profile), 
      shell=True
      )
