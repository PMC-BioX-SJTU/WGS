#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
from BasicWESorWGSLib import *
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Fastq Quality Control !!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='sample',action='store',help='sample info')
parser.add_argument('-o', required=True, dest='outpath4proj', action='store',help='outpath for the project')
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

def Define_QualityControlFq_Output():
	commands=['for(( i=0;i<${#fq1[@]};i++))\n',
		'do\n',
		'echo -e Time for quality control of fastq file was started at `date` "\\n" \\\n',
		'\t >>${outpath4step}/run.log\n',
		'gzip -dc ${inputdir}/${fq1[i]} |  \\\n',
		'\t${FASTX}/bin/fastx_quality_stats \\\n', 
		'\t-i -  ${index}\\\n',
		'\t-o ${outpath4step}/${name}.fq${i}.R1.stats &\n',
		'gzip -dc ${inputdir}/${fq2[i]} | \\\n',
		'\t${FASTX}/bin/fastx_quality_stats \\\n',
		'\t-i -  ${index}\\\n',
		'\t-o ${outpath4step}/${name}.fq${i}.R2.stats & \n',
		'wait\n\n',
		'python ${module_path}/Check_Read1_Read2.py \\\n',
		'\t-l ${outpath4step}/${name}.fq${i}.R1.stats \\\n',
		'\t-r ${outpath4step}/${name}.fq${i}.R2.stats\n\n'
		'Rscript ${module_path}/lineplot.R \\\n',
		'\t${outpath4step}/${name}.fq${i}.R1.stats \\\n',
		'\t${outpath4step}/${name}.fq${i}.Nu_distR1.png\n',

		'Rscript ${module_path}/lineplot.R \\\n',
		'\t${outpath4step}/${name}.fq${i}.R2.stats \\\n',
		'\t${outpath4step}/${name}.fq${i}.Nu_distR2.png \n\n' ,
		'${FASTX}/bin/fastq_quality_boxplot_graph.sh \\\n',
		'\t-i ${outpath4step}/${name}.fq${i}.R1.stats \\\n',
		'\t-o ${outpath4step}/${name}.fq${i}.QC_disR1.png \\\n',
		'\t-t ${name}.fq${i}.R1 &\n',
		'${FASTX}/bin/fastq_quality_boxplot_graph.sh \\\n',
		'\t-i ${outpath4step}/${name}.fq${i}.R2.stats \\\n',
		'\t-o ${outpath4step}/${name}.fq${i}.QC_disR2.png \\\n',
		'\t-t ${name}.fq${i}.R2 &\nwait\n',
		'done\n\n']
	return commands
	#output=sum([require_setting,log1,command,log2],[])
	#return (output)

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["FASTX"]
	param_opt=["index"]
	ref_opt=["ref"]
	in_opt=['fq1','fq2','inputtarget']
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='QualityControlFq'
        commands=Define_QualityControlFq_Output()
        
	### write shell scripts
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
