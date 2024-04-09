#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
from BasicWESorWGSLib import *
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Generating Germline Variants!')
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
	
def Define_GermlineVariantCalling_Output():
	commands=[
	'if [ ${datatype} == WGS ] ;then \n',
		'\tinterval=\' \'\n',
	'else\n',
		'\tinterval=\'-L ${bed} \'\n',
	'fi\n\n',
	'if [ ${vcftype} == gvcf ] ;then \n',
		'\tjava -jar ${GATK} -T HaplotypeCaller \\\n',
		'\t-R ${ref} -I ${inputtarget[0]} \\\n',
		'\t-o ${outpath4step}/${name}.raw.gvcf.gz -ERC GVCF \\\n',
		'\t-variant_index_type LINEAR -variant_index_parameter 128000 \n',
	'else\n',
		'\tjava -jar ${GATK} -T HaplotypeCaller \\\n',
		'\t-R ${ref} -I ${inputtarget} \\\n',
		'\t${interval}  -o ${outpath4step}/${name}.raw.vcf \\\n',
		'\t-stand_call_conf 30 --dbsnp ${dbsnp}  \n',
		#'\t-A RMSMappingQuality -A BaseCounts\n',
	'fi\n\n']
	return commands

if __name__ == '__main__':
	### requirements preparation
	soft_opt=["GATK"]
	param_opt=["vcftype","datatype"]
	ref_opt=["ref","bed","dbsnp"]
	in_opt=["inputtarget"]
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	#sec=['software','GermlineVariantCalling.parameter','reference','GermlineVariantCalling.output']
	step='GermlineVariantCalling'
        commands=Define_GermlineVariantCalling_Output()
	### write shell scripts
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
