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

	
def Define_DellyMergeSV_Output():
	command=[
		'${bcftools} merge   \\\n',
                '\t-m id -O b   \\\n',
                '\t-o ${outpath4step}/${name}.merged.INV.bcf   \\\n',
                '\t${inputtarget[@]} \n',

                '${bcftools} index     ${outpath4step}/${name}.merged.INV.bcf\n',

                '${delly} filter -t INV  \\\n',
                '-f germline  \\\n',
                '-o ${outpath4step}/${name}.germline.INV.bcf  \\\n',
                '${outpath4step}/${name}.merged.INV.bcf\n\n',

                '${bcftools} view  ${outpath4step}/${name}.germline.INV.bcf > ${outpath4step}/${name}.germline.INV.vcf\n\n'

                #'perl ${DellyFiltering} \\\n',
                #'\t-c ${Filtering_para} -n ${name} \\\n',
                #'\t-o ${outpath4step} \\\n',
                #'\t-g ${gaps} \\\n',
                #'\t-m ${centel} \n\n'

		]

	return (command)


if __name__ == '__main__':
	#requirements preparation
	soft_opt=["delly","bcftools"]
	param_opt=["sv_type"]
	ref_opt=["ref"]
	in_opt=["inputtarget"]
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='DellyMergeAllINV'
	commands=Define_DellyMergeSV_Output()
        ### write shell scripts
	Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
