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

	
def Define_CNVnatorCalling_Output():
	command=['export ROOTSYS=\"/home/zhujh/software/root\"\n',
		'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib\n\n',
		'${cnvnator}  -root ${outpath4step}/${name}.root \\\n',
		'\t-chrom ${chr} -unique \\\n',
		'\t-tree ${inputtarget[@]} \n\n',
		'geno_file=""\n',
		'for size in ${bin_range}\n',
		'do\n',
		'${cnvnator} -genome hg19 \\\n',
		'\t-root ${outpath4step}/${name}.root  \\\n',
		'\t-chrom ${chr} -his ${size} -d ${outpath4step}\n\n',

		'${cnvnator} -root ${outpath4step}/${name}.root  \\\n',
		'\t-chrom ${chr}    \\\n',
		'\t-stat ${size} \n\n',

		'${cnvnator} -root ${outpath4step}/${name}.root \\\n',
		'\t-chrom ${chr}  \\\n',
		'\t-partition ${size} -ngc \n',

		'${cnvnator} -root ${outpath4step}/${name}.root \\\n',
        	'\t-chrom ${chr} \\\n',
		'\t-call ${size} -ngc -unique >${outpath4step}/${name}.${size}.cnvnator \n\n',
		'geno_file="${geno_file} ${outpath4step}/${name}.${size}.cnvnator"\n\n',
		'done\n',
		'cat ${geno_file} >${outpath4step}/${name}.cnvnator\n',
		'perl ${CNVnatorFiltering}  \\\n',
		'\t-c ${Filtering_para}  \\\n', 
		'\t-n ${name} \\\n',
		'\t-g ${gaps}  \\\n',
		'\t-m ${centel} \\\n', 
		'\t-o ${outpath4step}\n\n']
	return (command)

def Makeclean_GermlineVariantCalling(path,name):
        clean=['clean:\n',
                '\trm'+path+name+'.raw*.vcf;\n',
                '\trm'+path+name+'.filtered*.vcf;\n']
        retrun (clean)

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["cnvnator"]
	param_opt=["Filtering_para","bin_range","chr"]
	ref_opt=["ref","gaps","centel"]
	in_opt=["inputtarget"]
	out_opt=["outtarget"]
	inhouse_opt=["module_path","CNVnatorFiltering"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='CNVnatorCalling'
	commands=Define_CNVnatorCalling_Output()
	### write shell scripts
	Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)

