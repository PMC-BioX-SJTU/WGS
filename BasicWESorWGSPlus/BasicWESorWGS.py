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
parser.add_argument('-s', required=True, dest='sample',action='store',help='sample info')
parser.add_argument('-m', required=True, dest='outfile', action='store',help='outfile for make')
parser.add_argument('-o', required=True, dest='outpath', action='store',help='outpath for project')

if len(sys.argv) <= 2: 
	parser.print_help() 
	sys.exit(1) 
else: 
	args = parser.parse_args()
config = args.config
outfile = args.outfile

sample = args.sample
outpath=args.outpath
from BasicWESorWGSLib import *

def Read_config_info(config):
	actSteps=[]
	ConfigInfo=Read_ConfigFile(config)
	module_path=ConfigInfo.Get_OptionItems('inhouse_scripts','module_path')
	steps=ConfigInfo.Get_OptionItems('actSteps','StepstoRun')
	actSteps=steps.split(',')
	return(module_path,actSteps)

def Ext_modules(sample,module_path,actSteps,outpath,outfile):
	SampleInfo=Obtain_SampleInfo(sample)
	SampleInfo.Check_Length()
	SampleInfo.Define_Names()
	SampleInfo.Define_Dirs()
	SampleInfo.Define_FQs()
	SamInfo=SampleInfo
	for j in range(0,len(actSteps)):
		tar=[]
		ext_out=[]
		for i in range(0,len(SamInfo.names)):
			### prepare input info
			name=SamInfo.names[i]
			fq=SamInfo.fqs[i]
			directory=SamInfo.dirs[i]
			if(actSteps[j]== 'QualityControlFq' or actSteps[j]== 'MappingAndCleanUp' ):
				sam_info='\''+name+ '|' + fq +':'+directory +'\''
			else:
				sam_info='\''+name+'\''
			script=module_path+ '/' + actSteps[j] + '.py'
			### commands act steps
			cmd='python '+ script + ' -c '+ config + ' -s ' + sam_info + ' -o ' + outpath + ' -m BasicWESorWGS'
			os.system(cmd)
			### command for make file
			log='$(outpath)/'+name+'/'+actSteps[j]+'/'+actSteps[j]+'.log'
			tar.append(log)
			ext_out.append('$(eval $(call steps,'+name+ ','+log+'))\n')
		tar_out=' \\\n'.join(tar)
		Write_makefile(outfile,outpath+'/BasicWESorWGS',tar_out,ext_out,actSteps[j])

if __name__ == '__main__':
	(module_path,actSteps)= Read_config_info(config)
	Ext_modules(sample,module_path,actSteps,outpath,outfile)


