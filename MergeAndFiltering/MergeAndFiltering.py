#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
sys.path.append('/home/zhujh/Codes/Pipelines/BasicWESorWGSPlus')
sys.path.append('/home/zhujh/Codes/Pipelines/FamilyBasic')
from BasicWESorWGSLib import *
from FamilyBasicLib import *
import re
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='ped',action='store',help='ped info')
parser.add_argument('-o', required=True, dest='outpath', action='store',help='outpath for project')
parser.add_argument('-m', required=True, dest='outfile', action='store',help='outfile for make')
if len(sys.argv) <= 2: 
	parser.print_help() 
	sys.exit(1) 
else: 
	args = parser.parse_args()
config = args.config
sample = args.ped
outpath=args.outpath
outfile = args.outfile


def Read_config_info(config):
	actSteps=[]
	cf.read(config)
	module_path=cf.get('inhouse_scripts','module_path')
	actSteps=cf.get('actSteps','StepstoRun')
	waitSteps=cf.get('actSteps','StepstoWait')
	return(module_path,actSteps,waitSteps)

def Ext_modules(ped,module_path,actSteps,waitSteps,outpath,outfile,chrGroup):
	FamilyInfo=Obtain_FamilyInfo(ped)
	FamilyInfo.Get_famGroup(30)
	FamInfo=FamilyInfo
	combine_tar=[]
	combine_extOut=[]
	genotype_tar=[]
	genotype_extOut=[]
	vcf=[]
	for j in range(0, len(chrGroup)):
		gvcf=[]
		for i in range(0,len(FamInfo.famGroup)):
			member=':'.join(FamInfo.famGroup[i])
			script=module_path+ '/CombineGVCFs.py'
			name='FamGroup'+str(i)+'_'+chrGroup[j]
			sample=name+'|'+member
			cmd1='python ' + script + ' -c '+ config + ' -s ' + '\''+ sample +'\'' + ' -o ' + outpath + ' -m MergeAndFiltering'
			os.system(cmd1)
			### make input
			log='$(outpath)/' + name +'/CombineGVCFs/CombineGVCFs.log'
			#gvcf.append('${InputPath}/'+name+'/CombineGVCFs/'+name+'.raw.'+chrGroup[j]+'.gvcf.gz')
			gvcf.append(name)
			combine_tar.append(log)
			combine_extOut.append('$(eval $(call steps,'+name+ ','+log+'))\n')
		gvcf_lst=':'.join(gvcf)
		gvcf_sam=chrGroup[j]+'|'+gvcf_lst
		script=module_path+ '/GenotypeGVCFs.py'
		cmd1='python ' + script + ' -c '+ config + ' -s ' + '\''+ gvcf_sam +'\'' + ' -o ' + outpath + ' -m MergeAndFiltering'
		os.system(cmd1)	
		log1='$(outpath)/' + chrGroup[j] +'/GenotypeGVCFs/GenotypeGVCFs.log'
		#vcf.append('${InputPath}/'+chrGroup[j]+'/GenotypeGVCFs/'+chrGroup[j]+'.raw.vcf.gz')
		genotype_tar.append(log1)
		genotype_extOut.append('$(eval $(call steps,'+chrGroup[j]+ ','+log1+'))\n')
	vcf_lst=':'.join(chrGroup)
	vcf_all='all|'+vcf_lst
	script=module_path+ '/CombineVariants.py'
	cmd1='python ' + script + ' -c '+ config + ' -s ' + '\''+ vcf_all +'\'' + ' -o ' + outpath + ' -m MergeAndFiltering'
	os.system(cmd1)
	### command for make file	
	combine_tarOut=' \\\n'.join(combine_tar)
	Write_makefile(outfile,outpath+'/MergeAndFiltering',combine_tarOut,combine_extOut,'CombineGVCFs')
	genotype_tarOut=' \\\n'.join(genotype_tar)	
	Write_makefile(outfile,outpath+'/MergeAndFiltering',genotype_tarOut,genotype_extOut,'GenotypeGVCFs')
	### merge QCs
	all_mem=[]
	for i in range(0,len(FamInfo.famGroup)):
		all_mem.append(':'.join(FamInfo.famGroup[i]))	
	members='all|'+':'.join(all_mem)
	#script=module_path+ '/MergeGVCFsOther.py'
	#cmd1='python ' + script + ' -c '+ config + ' -s ' + '\''+ members +'\'' + ' -o ' + outpath + ' -m MergeAndFiltering'
	#os.system(cmd1)

	script=module_path+ '/MergeQCs.py'
	cmd1='python ' + script + ' -c '+ config + ' -s ' + '\''+ members +'\'' + ' -o ' + outpath + ' -m MergeAndFiltering'	
	os.system(cmd1)
	script=module_path+ '/VariantFiltering.py'
	cmd1='python ' + script + ' -c '+ config + ' -s  all  -o ' + outpath + ' -m MergeAndFiltering'
	os.system(cmd1)
	script=module_path+ '/QualityControlVCF.py'
	cmd1='python ' + script + ' -c '+ config + ' -s  all  -o ' + outpath + ' -m MergeAndFiltering'
        os.system(cmd1)	


if __name__ == '__main__':
	chrGroup=['chr1','chr2','chr3','chr4',
		'chr5','chr6','chr7','chr8','chr9', 
		'chr10','chr11','chr12','chr13','chr14',
		'chr15','chr16','chr17','chr18','chr19',
		'chr20','chr21','chr22','chrX','chrY','chrM']
	(module_path,actSteps,waitSteps)= Read_config_info(config)
	Ext_modules(sample,module_path,actSteps,waitSteps,outpath,outfile,chrGroup)
