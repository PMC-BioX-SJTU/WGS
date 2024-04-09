#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
import re
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='ped',action='store',help='ped info')
parser.add_argument('-m', required=True, dest='outfile', action='store',help='outfile for make')
parser.add_argument('-o', required=True, dest='outpath', action='store',help='outpath for project')

if len(sys.argv) <= 2: 
	parser.print_help() 
	sys.exit(1) 
else: 
	args = parser.parse_args()
config = args.config
outfile = args.outfile
sample = args.ped
outpath=args.outpath

sys.path.append('/home/zhujh/Codes/Pipelines/BasicWESorWGSPlus')
from BasicWESorWGSLib import *
sys.path.append('/home/zhujh/Codes/Pipelines/FamilyBasic')
from FamilyBasicLib import *

def Read_config_info(config):
	actSteps=[]
	ConfigInfo=Read_ConfigFile(config)
	module_path=ConfigInfo.Get_OptionItems('inhouse_scripts','module_path')
	steps=ConfigInfo.Get_OptionItems('actSteps','StepstoRun')
	actSteps=steps.split(',')
	return(module_path,actSteps)

def Ext_modules(ped,module_path,actSteps,outpath,outfile):
	FamilyInfo=Obtain_FamilyInfo(ped)
	FamilyInfo.Get_famMem()
	FamilyInfo.Output_famPed(outpath+'/FamilyBasic')
	FamInfo=FamilyInfo
	combine_tar=[]
	combine_extOut=[]
	combine_tar1=[]
	combine_extOut1=[]
	for fam in FamInfo.famMem:
		## comphet
		script=module_path+ '/CompHetFinderAndFiltering.py'	
		cmd1='python ' + script + ' -c ' + config + ' -s '+ '\''+fam+"|"+'\'' + ' -o ' + outpath + ' -m FamilyBasic'
		os.system(cmd1)
		log='$(outpath)/' + fam +'/CompHetFinderAndFiltering/CompHetFinderAndFiltering.log'
		combine_tar.append(log)
		combine_extOut.append('$(eval $(call steps,'+fam+ ','+log+'))\n')
		## ADARXlink
		script=module_path+ '/ADARXlinkFinderAndFiltering.py'
		cmd1='python ' + script + ' -c ' + config + ' -s '+ '\''+fam+"|"+'\'' + ' -o ' + outpath + ' -m FamilyBasic'
		os.system(cmd1)
		log1='$(outpath)/' + fam +'/ADARXlinkFinderAndFiltering/CompHetFinderAndFiltering.log'
		combine_tar1.append(log1)
		combine_extOut1.append('$(eval $(call steps,'+fam+ ','+log1+'))\n')
	## comphet make
	combine_tarOut=' \\\n'.join(combine_tar)
	Write_makefile(outfile,outpath+'/FamilyBasic',combine_tarOut,combine_extOut,'CompHetFinderAndFiltering')
	## ADARXlink make
	combine_tarOut1=' \\\n'.join(combine_tar1)
	Write_makefile(outfile,outpath+'/FamilyBasic',combine_tarOut,combine_extOut,'ADARXlinkFinderAndFiltering')	
	##DeNovo mutation calling
	script=module_path+ '/DeNovoFinderAndFiltering.py'
	cmd1='python ' + script + ' -c ' + config + ' -s '+ '\'all|\'' + ' -o ' + outpath + ' -m FamilyBasic'
	os.system(cmd1)
	
if __name__ == '__main__':
	(module_path,actSteps)= Read_config_info(config)
	Ext_modules(sample,module_path,actSteps,outpath,outfile)


