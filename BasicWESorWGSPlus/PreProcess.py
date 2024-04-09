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
parser.add_argument('-s', required=True, dest='sample',action='store',help='sample info')
parser.add_argument('-o', required=True, dest='outpath4proj', action='store',help='outpath for the project')

if len(sys.argv) <= 3: 
	parser.print_help() 
	sys.exit(1) 
else: 
	args = parser.parse_args()
config = args.config
outpath4proj = args.outpath4proj
sample = args.sample

def Input_SampleAndOutpath(tar,sample,step,outpath4proj):
	if (re.findall('.gz',sample)):
		(name,fq)=sample.split('|')
		### mkdir for this step
		outpath4step=outpath4proj+'/'+name+'/'+step
		if  not os.path.exists(outpath4step):
			os.makedirs(outpath4step)
		### determine input target file
		if tar != 'NO':
			samp_out=[name,tar,outpath4proj,outpath4step]
		else:
			samp_out=[name,fq,outpath4proj,outpath4step]
	else:
		name=sample
		outpath4step=outpath4proj+'/'+step
		samp_out=[name,tar,outpath4proj,outpath4step]
		if  not os.path.exists(outpath4step):
			os.makedirs(outpath4step)
	return(samp_out)

def Input_Requiremens(config,sec,opt,sample,step,outpath4proj):
	cf.read(config)
	requirements=[]
	tar=cf.get(step + ".input","input.target")
	requirements.append(Input_SampleAndOutpath(tar,sample,step,outpath4proj))
	for i in range(0,len(sec)):
		setting=[]
		for j in range(0,len(opt[i])):
        		setting.append(cf.get(sec[i], opt[i][j]))
		requirements.append(setting)
        return(requirements)

def Check_SampleRawFq(name,fq):
	(r1,r2,inputdir)=fq.split(":")
	raw_r1=r1.split()
	raw_r2=r2.split()
	if(len(raw_r1) != len(raw_r2)):
		print name,"left Read cannot match right Read completely!!"
		sys.exit(1)
	else:
		return(raw_r1,raw_r2,r1,r2,inputdir)

def Check_Exists(file):
        for i in range(0,len(file)):
                if not os.path.exists(file[i]):
                        print "We cannot find the file "+file[i]+"!!"
                        sys.exit(1)

def Define_PreProcess_Output(requirements):
	(raw_r1,raw_r2,r1,r2,inputdir)=Check_SampleRawFq(requirements[0][0],requirements[0][1])
	if(len (raw_r1) ==1 ):
		out=['echo -e The Step PreProcess started at `date` "\\n" >>'+ requirements[0][3]+'/run.log\n',
		     'mv ' +inputdir+'/'+ r1 + ' ' + requirements[0][3]+'/' + requirements[0][0] + '.R1.fq.gz\n', 
		     'mv ' +inputdir+'/'+ r2 + ' ' + requirements[0][3]+'/' + requirements[0][0] + '.R2.fq.gz\n',
		     'echo -e The Step PreProcess ended at `date` "\\n" >>'+ requirements[0][3]+'/run.log\n',
		     'mv '+requirements[0][3]+'/run.log'+' '+ requirements[0][3]+'/PreProcess.log;\n']
	else:
		#for i in range(0,len(raw_r1)):
		#	raw_r1[i]=inputdir+'/'+raw_r1[i]
		#	raw_r2[i]=inputdir+'/'+raw_r2[i]
		out=['echo -e The Step PreProcess started at `date` "\\n" >>'+ requirements[0][3]+'/run.log\n',
		    'cat ' + ' '.join([inputdir+'/'+x for x in raw_r1]) + '|' + ' gzip -d '+'|' + ' gzip - >' + requirements[0][3]+'/' + requirements[0][0] + '.R1.fq.gz\n',
		    'cat ' + ' '.join([inputdir+'/'+x for x in raw_r2]) + '|' + ' gzip -d '+'|' + ' gzip - >' + requirements[0][3]+'/' + requirements[0][0] + '.R2.fq.gz\n',
		     'echo -e The Step PreProcess ended at `date` "\\n" >>'+ requirements[0][3]+'/run.log\n',
		     'mv '+requirements[0][3]+'/run.log'+' '+ requirements[0][3]+'/PreProcess.log\n']

	output=sum([out],[])
	return (output)

def Makeclean_MappingAndCleanUp(path,name):
	clean=['clean:\n',
		'\trm'+path+name+'.align.*;\n',
		'\trm'+path+name+'.realigneri.*;\n',
		'\trm'+path+name+'*recal.table\n']
	retrun (clean)

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["module_path"]
	param_opt=["merge"]
	ref_opt=["ref"]
	opt=[soft_opt,param_opt,ref_opt]
	sec=['software','PreProcess.parameter','reference']
	### SampleAndOutpat, software, parameter,reference
	requirements=Input_Requiremens(config,sec,opt,sample,"PreProcess",outpath4proj)
	
	### write makefile
	makefilename=requirements[0][3]+'/PreProcess.sh'
	f=open(makefilename,"w")
	output=Define_PreProcess_Output(requirements)
	f.writelines(output)
	f.close
	
