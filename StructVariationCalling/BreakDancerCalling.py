#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
sys.path.append('/home/zhujh/Codes/Pipelines/BasicWESorWGSPlus')
from BasicWESorWGSLib import *
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='sample',action='store',help='ped file info')
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

	
def Define_BreakDancerCalling_Output():
	### setting paramenter
        ### Merging sample lists of vcfs or bams
	
	breakdancer=['perl ${bam2cfg} \\\n',
		'\t${inputtarget[@]} > ${outpath4step}/all.list\n',
		'${breakdancer}  -m 1000000   \\\n',
		'-h  ${outpath4step}/all.list >${outpath4step}/${name}.SV.output\n\n',
		'egrep \'DEL|INV|INS\' ${outpath4step}/${name}.SV.output >${outpath4step}/${name}.DEL.INV.INS.output\n',
		'perl ${BreakDancerFiltering}  \\\n',
		'\t-c  ${Filtering_para}  \\\n',
		'\t-n ${name}  -g ${gaps} \\\n',
		'\t-m ${centel} \\\n',
		'\t-o ${outpath4step}\n\n',
	]


	return (breakdancer)

def Makeclean_GermlineVariantCalling(path,name):
        clean=['clean:\n',
                '\trm'+path+name+'.raw*.vcf;\n',
                '\trm'+path+name+'.filtered*.vcf;\n']
        retrun (clean)

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["bam2cfg","breakdancer"]
	param_opt=["Filtering_para"]
	ref_opt=["ref","gaps","centel"]
	in_opt=['inputtarget']
	out_opt=["outtarget"]
	inhouse_opt=["module_path","BreakDancerFiltering"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
        step='BreakDancerCalling'
        commands=Define_BreakDancerCalling_Output()
        ### write shell scripts
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)

