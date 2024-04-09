#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
sys.path.append('/home/zhujh/Codes/Pipelines/BasicWESorWGSPlus/')
from BasicWESorWGSLib import *
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='sample',action='store',help='sample info')
parser.add_argument('-o', required=True, dest='outpath4proj', action='store',help='outfile for make')
parser.add_argument('-m', required=True, dest='module', action='store',help='module name')

    
if len(sys.argv) <= 3: 
	parser.print_help() 
	sys.exit(1) 
else: 
	args = parser.parse_args()
config = args.config
outpath4proj = args.outpath4proj
sample = args.sample
module=args.module

	
def Define_DellyCalling_Output():
	command=['for sv in ${sv_type}\n',
		'do\n'
		'\t${delly} call -t ${sv} \\\n',
		'\t\t-g ${ref} \\\n',
		'\t\t-o ${outpath4step}/${name}.${sv}.bcf \\\n',
		'\t\t-x ${hg19_excl} ${inputtarget[@]}\n',
		'done\n']		

	return (command)


if __name__ == '__main__':
	#requirements preparation
	soft_opt=["delly","bcftools"]
	param_opt=["Filtering_para","sv_type"]
	ref_opt=["ref",'hg19_excl',"gaps","centel"]
	in_opt=["inputtarget"]
	out_opt=["outtarget"]
	inhouse_opt=["module_path","DellyFiltering"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='DellyCalling'
	commands=Define_DellyCalling_Output()
        ### write shell scripts
	Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
