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

	
def Define_PindelCalling_Output():
	command=['for(( i=0;i<${#inputtarget[@]};i++))\n',
		'do\n',
		'bam=`basename ${inputtarget[i]}`\n',
		'bam=${bam%%.*}\n',
		'perl ${module_path}/extract_insertsize.pl \\\n',
		'\t${name} ${InputPath}/${bam}/QualityControlBam/${bam}.insert_size.txt \\\n',
		'\t${inputtarget[i]} >${outpath4step}/${bam}.config\n',
		'${pindel} -w 0.1 -x 5 -B 0 -T 4 \\\n',
		'\t-f ${ref} \\\n', 
		'\t-i ${outpath4step}/${bam}.config \\\n',
		'\t-c ALL -o ${outpath4step}/${bam} \n\n',
		'done\n',
		'for type in ${sv_type}\n',
		'do\n',
		'geno_file="" \n',
		'for(( i=0;i<${#inputtarget[@]};i++))\n',
		'do\n',
		'bam=`basename ${inputtarget[i]}`\n',
                'bam=${bam%%.*}\n',
		'${pindel2vcf} -p ${outpath4step}/${bam}_${type} \\\n',
		'\t-r ${ref} \\\n',
		'\t-R UCSC_hg19  \\\n',
		'\t-d 20101123 \\\n',
		'\t-v ${outpath4step}/${bam}_${type}.vcf -G \n\n',
		'${bgzip} ${outpath4step}/${bam}_${type}.vcf \n',
		'${tabix} -p vcf ${outpath4step}/${bam}_${type}.vcf.gz \n',
		'geno_file="${geno_file} ${outpath4step}/${bam}_${type}.vcf.gz"\n',
		'done\n',
		'${bcftools} merge  \\\n',
		'\t-m id -O v  \\\n',
		'\t-o ${outpath4step}/${name}.merged.${type}.vcf  \\\n',
		'\t${geno_file} \n',
		'done\n',
		'perl ${Module_path}/PindelFiltering.pl  \\\n', 
		'\t-c ${Module_path}/StructVariationFitering.config  \\\n',
		'\t-g ${gaps} \\\n',
		'\t-m ${centel}  \\\n',
		'\t-n ${name}  \\\n',
		'\t-o ${outpath4step}\n\n']

	#output=sum([require_setting,log1,command,log2],[])
	return (command)

def Makeclean_GermlineVariantCalling(path,name):
        clean=['clean:\n',
                '\trm'+path+name+'.raw*.vcf;\n',
                '\trm'+path+name+'.filtered*.vcf;\n']
        retrun (clean)

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["pindel","bgzip","tabix"]
	param_opt=["Filtering_para","sv_type"]
	ref_opt=["ref","gaps","centel"]
	in_opt=["inputtarget"]
        out_opt=["outtarget"]
	inhouse_opt=["module_path","PindelFiltering"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='PindelCalling'
	commands=Define_PindelCalling_Output()
        ### write shell scripts
	Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
