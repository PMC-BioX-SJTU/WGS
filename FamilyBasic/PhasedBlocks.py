#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
sys.path.append('/opt/NfsDir/UserDir/sunfl/Codes/Pipelines/Modules/')
from PreProcess import *
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='sample',action='store',help='sample info')
parser.add_argument('-o', required=True, dest='outpath4proj', action='store',help='outfile for make')


    
if len(sys.argv) <= 3: 
	parser.print_help() 
	sys.exit(1) 
else: 
	args = parser.parse_args()
config = args.config
outpath4proj = args.outpath4proj
sample = args.sample
  

	
def Define_VariantFiltering_Output(requirements):
	sample=['name=',requirements[0][0],'\n','outpath4proj=',requirements[0][2],'\n','outpath4step=',requirements[0][3],'\n','inputtarget=',requirements[0][1],'\n\n']
	software=['module_path=',requirements[1][0],'\n\n']
	parameter=['number_of_nonfounders=',requirements[2][0],'\n','inclusion_order=',requirements[2][1],'\n\n']
	reference=['gap_file=',requirements[3][0],'\n\n']

	log=['echo -e The Step VariantFiltering started at `date` "\\n" >>${outpath4step}/run.log;\n',
		'for file in ${inputtarget}\n',
                'do\n',
                'if [ ! -f "${file}" ]; then\n',
                '\techo -e the target file ${file} cannot bed generated successfully!!\n'
                '\texit 0\n',
                'fi\n',
                'done\n']

	vcf2csv=['perl vcf2ISBcsv.pl ${inputtarget} >${outpath4step}/all.final.csv']
	hmm=['perl ${module_path}/suttonian_hmm.pl \\\n',
		'\t--number_of_nonfounders=${number_of_nonfounders} \\\n',
		'\t--file_to_analyze=${outpath4step}/all.final.csv \\\n', 
		'\t--inclusion_order=${inclusion_order} \\\n',
		'\t--output_dir=${outpath4step} \\\n',
		'\t--project=${name} \\\n',
		'\t--pedigree_name=${name} \\\n',
		'\t--exclude_fully_heterozygous_vectors \\\n',
		'\t--exclude_partially_called_vectors\n\n']

	block_nibbler=['perl  ${module_path}/suttonian_block_nibbler.pl  \\\n',
		'\t--infile_directory=${outpath4step} \n\n']
	statistics=['perl  ${module_path}/statistics_for_HMM_Suttonian_blocks.pl  \\\n',
		'\t--reference_genome=hg19 \\\n', 
		'\t--infile_directory=${outpath4step} \\\n',
		'\t--gap_file=${gap_file} \n\n']
	smooth=['perl ${module_path}/parse_and_smooth_HMM_Suttonian_blocks.pl \\\n',
		'\t--infile_directory=${outpath4step}\n']

	log1=['echo -e The Step VariantFiltering ended at `date` "\\n" >>${outpath4step}/run.log;\n',
		'for file in ${outpath4step}/${outtarget}\n',
                'do\n',
                'if [ ! -f "${file}" ]; then\n',
                '\techo -e the target file ${file} cannot bed generated successfully!!\n'
                '\texit 0\n',
                'fi\n',
                'done\n']

	output=sum([sample,software,reference,parameter,log,vcf2csv,hmm,block_nibbler,statistics,smooth,log1],[])
	return (output)


if __name__ == '__main__':
	#requirements preparation
	soft_opt=["module_path"]
	param_opt=["number_of_nonfounders","inclusion_order"]
	ref_opt=["gap_file"]
	out_opt=["output.target"]
	opt=[soft_opt,param_opt,ref_opt,out_opt]
	sec=['software','PhasedBlocks.parameter','reference','PhasedBlocks.output']
	### SampleAndOutpat, software, parameter,reference
	requirements=Input_Requiremens(config,sec,opt,sample,"PhasedBlocks",outpath4proj)
	### write makefile
	makefilename=requirements[0][3]+'/PhasedBlocks.sh'
	f=open(makefilename,"w")
	output=Define_VariantFiltering_Output(requirements)
	f.writelines(output)
	f.close
