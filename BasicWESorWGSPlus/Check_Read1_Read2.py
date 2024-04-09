#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for checking count of read1 and read2!')
parser.add_argument('-l', required=True, dest='stat1',action='store',help='config file for the module')
parser.add_argument('-r', required=True, dest='stat2',action='store',help='sample info')

if len(sys.argv) <= 2:
        parser.print_help()
        sys.exit(1)
else:
        args = parser.parse_args()
stat1 = args.stat1
stat2 = args.stat2

def Check_Read1_Read2(stat1,stat2):
        f1=open(stat1,'r')
        lines1=f1.readlines()
        items1=lines1[1].split()
        f1.close()
        f2=open(stat2,'r')
        lines2=f2.readlines()
        items2=lines2[1].split()
        f2.close()
        if(items1[1] != items2[1]):
                print "Read1 count is not same ad read2 count! "
                sys.exit(1)

if __name__ == '__main__':
	Check_Read1_Read2(stat1,stat2)
