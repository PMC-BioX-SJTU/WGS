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

def Ext_modules(ped,module_path,actSteps,outpath,outfile):
	FamilyInfo=Obtain_FamilyInfo(ped)
	FamilyInfo.Get_famMem()
	FamInfo=FamilyInfo
	for j in range(0,len(actSteps)):
		tar=[]
		ext_out=[]
		famList=[]
		i=0
		for fam in FamInfo.famMem:
			famList.append(fam)
			### prepare input info
			bam=':'.join(FamInfo.famMem[fam])
			sam_info='\'' + fam+ '|' + bam  +'\''
			script=module_path+ '/' + actSteps[j] + '.py'
			### commands act steps
			cmd='python '+ script + ' -c '+ config + ' -s ' + sam_info + ' -o ' + outpath + ' -m StructVariationCalling'
			os.system(cmd)
			### command for make file
			log='$(outpath)/StructVariationCalling/' + fam +'/'+actSteps[j]+'/' + actSteps[j]+'.log'
			tar.append(log)
			ext_out.append('$(eval $(call steps,'+fam+ ','+log+'))\n')
			i+=1
		tar_out=' \\\n'.join(tar)
		Write_makefile(outfile,outpath+'/StructVariationCalling',tar_out,ext_out,actSteps[j])
	### DellyMergeSV
	#### merge DEL site
	famListStr=':'.join(famList)
	sam_info='\'all|' + famListStr  +'\''
	script=module_path+ '/DellyMergeDEL.py'
	cmd='python '+ script + ' -c '+ config + ' -s ' + sam_info + ' -o ' + outpath + ' -m StructVariationCalling'
	os.system(cmd)
	#### merge DEL bcf
	script=module_path+ '/DellyMergeAllDEL.py'
        cmd='python '+ script + ' -c '+ config + ' -s ' + sam_info + ' -o ' + outpath + ' -m StructVariationCalling'
        os.system(cmd)
	#### DUP
	script=module_path+ '/DellyMergeDUP.py'
        cmd='python '+ script + ' -c '+ config + ' -s ' + sam_info + ' -o ' + outpath + ' -m StructVariationCalling'
        os.system(cmd)
	#### merge DUP bcf
	script=module_path+ '/DellyMergeAllDUP.py'
	cmd='python '+ script + ' -c '+ config + ' -s ' + sam_info + ' -o ' + outpath + ' -m StructVariationCalling'
	os.system(cmd)
	#### INV
	script=module_path+ '/DellyMergeINV.py'
        cmd='python '+ script + ' -c '+ config + ' -s ' + sam_info + ' -o ' + outpath + ' -m StructVariationCalling'
        os.system(cmd)
	#### merge INV bcf
	script=module_path + '/DellyMergeAllINV.py'
	cmd='python '+ script + ' -c '+ config + ' -s ' + sam_info + ' -o ' + outpath + ' -m StructVariationCalling'
	os.system(cmd)

if __name__ == '__main__':
	(module_path,actSteps)= Read_config_info(config)
	Ext_modules(sample,module_path,actSteps,outpath,outfile)


