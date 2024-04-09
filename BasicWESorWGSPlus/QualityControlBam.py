#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
from BasicWESorWGSLib import *
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
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


def Define_QualityControlBam_Output():
	cross_contam=['${verifyBamID} --vcf ${EAS_1000g} \\\n', 
		'--bam ${inputtarget}  \\\n',
		'--out ${outpath4step}/${name}.verifybamid  \\\n',
		'--verbose -ignoreRG\n\n']
	hs=['if [ ${datatype} == WGS ] ;then\n',
		'java -jar ${picard} CollectRawWgsMetrics\\\n', 
		'\tI=${inputtarget}  \\\n',
		'\tO=${outpath4step}/${name}.hs.raw  \\\n',
		'\tR=${ref} INCLUDE_BQ_HISTOGRAM=true\\\n',
		'\tCOVERAGE_CAP=250 MINIMUM_BASE_QUALITY=20 MINIMUM_MAPPING_QUALITY=20\n\n',
		'else\n',
		'java -jar ${picard}  CollectHsMetrics \\\n',
		'\tI=${inputtarget}  \\\n',
		'\tO=${outpath4step}/${name}.hs.raw  \\\n',
		'\tR=${ref} BAIT_INTERVALS=${bed} \\\n',
		'\tTARGET_INTERVALS=${bed}\n\n',
		'fi\n\n']

	alignment=['java -jar ${picard}  CollectAlignmentSummaryMetrics \\\n \
	R=${ref} I=${inputtarget}  \\\n \
	O=${outpath4step}/${name}.alignment.raw\n\n']

	qy=['java -jar ${picard}  CollectQualityYieldMetrics \\\n \
	R=${ref} I=${inputtarget} \\\n \
	O=${outpath4step}/${name}.qualigyYield.raw\n\n']

	insert_size=['java -jar ${picard}  CollectInsertSizeMetrics \\\n \
	R=${ref} I=${inputtarget}  \\\n \
	O=${outpath4step}/${name}.insert_size.txt  \\\n \
	H=${outpath4step}/${name}.insert_size.pdf M=0.5 \n',
	'convert  ${outpath4step}/${name}.insert_size.pdf \\\n \
	${outpath4step}/${name}.insert_size.png\n\n']
	
	qc_bam=['perl ${module_path}/QCBam_table.pl ${name} \\\n \
	${outpath4step}/${name}.verifybamid.selfSM \\\n \
	${outpath4step}/${name}.hs.raw \\\n \
	${outpath4step}/${name}.alignment.raw \\\n \
	${outpath4step}/${name}.qualigyYield.raw ${datatype}\n',
	'perl ${module_path}/QCBam_plot.pl \\\n \
	${outpath4step}/${name}.hs.raw  \\\n \
	${name}.depth.freq ${name}.cumu.freq ${name}\n',
	'rm  ${outpath4step}/${name}.verifybamid.depthRG ${outpath4step}/${name}.verifybamid.log ${outpath4step}/${name}.verifybamid.selfRG ${outpath4step}/${name}.verifybamid.depthSM ${outpath4step}/${name}*freq ${outpath4step}/${name}*R\n\n']	
	commands=sum([cross_contam, hs, alignment, qy, insert_size,qc_bam],[])
	return commands

if __name__ == '__main__':
	#requirements preparation
	soft_opt=[ "verifyBamID","picard"]
	param_opt=["datatype"]
	ref_opt=["ref","bed","EAS_1000g"]
	in_opt=['inputtarget']
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='QualityControlBam'
        commands=Define_QualityControlBam_Output()
        ### write shell scripts
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
