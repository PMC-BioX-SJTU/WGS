#/usr/bin/python
import ConfigParser
import string
import os
import argparse
import sys
import re
cf = ConfigParser.ConfigParser()
class Read_ConfigFile:
	def __init__(self,FileName):
		cf.read(FileName)
	
	def Get_Options(self,sec_name):
		info=cf.options(sec_name)
		return info

	def Get_OptionItems(self,sec_name,opt_name):
		item=cf.get(sec_name, opt_name)
		return item
	
	def Get_OptionLists(self,sec,opt):
		lists=[]
		for k in range(0,len(opt)):
			lists.append(cf.get(sec, opt[k]))
		return lists
	
	def Get_OptionLists4shell(self,sec,opt,i):
		setting=['###'+sec[i]+'\n']
		for k in range(0,len(opt[i])):
			setting.append(opt[i][k]+'='+cf.get(sec[i], opt[i][k])+'\n')
		setting.append('\n')
		return setting


class Obtain_SampleInfo(Read_ConfigFile):
	def __init__(self,sample): 
		self.sample=sample
		Read_ConfigFile.__init__(self,self.sample)
	
	def Get_Info(self,sec_name):
		info=self.Get_Options(sec_name)
		Length=len(info)
		return(info,Length)
        
	def Check_Length(self):
		(self.SampleName,len1)=self.Get_Info("Sample.name")
		(self.SampleData,len2)=self.Get_Info("Sample.data")
		(self.SampleDir, len3)=self.Get_Info("Sample.dir")
		if (len1 != len2 or len1 != len3 ): ### check sample.ini input
        	        print "Sample count is not matched with fq count!!!"
                	sys.exit(1)

        def Define_Names(self):
                self.names=self.Get_OptionLists("Sample.name",self.SampleName)
	
	def Define_Dirs(self):
		self.dirs=self.Get_OptionLists("Sample.dir",self.SampleDir)

	def Define_FQs(self):
		self.fqs=self.Get_OptionLists("Sample.data",self.SampleData)

class InOutput_Target_Pattern():
	def __init__(self,index):
		self.index=index
	
	def Split_Str(self,string,sep):
		lst=string.split(sep)
		return lst
	
	def Target_Simp(self,Path_lst,Mem,sep):
		if(len(Path_lst) != len(Mem)):
			print name,"Path count cannot match Memebers completely!!"
			sys.exit(1)
		tar_lst=[]
		for k in range(0,len(Path_lst)):
			Lst=self.Split_Str(Mem[k],sep)
			Target=[Path_lst[k]+ x  for x in Lst]
			tar_lst.append(Target)
		return(tar_lst)

	def Target_Lst(self,InputTarget,Mem_lst,sep):
		lst=InputTarget.split("|")
 		suf_lst=lst[2].split(":")
		Path_lst=[lst[0]]*len(suf_lst)
		step_lst=[lst[1]]*len(suf_lst)
 		Mem=[Mem_lst]*len(suf_lst)
		len1=len(Mem)
		if(len(Path_lst) != len1 or  len(step_lst) != len1 or  len(suf_lst) != len1):
			print name,"some list count cannot match Memebers completely!!"
			sys.exit(1)
		tar_lst=[]
		for j in range(0,len(Path_lst)):
			list1=self.Split_Str(Mem[j],sep)
			if(len(lst) == 3):
				Target=[Path_lst[j]+ x + step_lst[j] + x + suf_lst[j]  for x in list1]
			else:
				Target=[Path_lst[k]+ x  for x in Lst]
			tar_lst.append(Target)
		return(tar_lst)			
	
	def Target_Join(self,sep_lst,tar_lst,head_lst):
		if(len(sep_lst) != len(tar_lst) or len(sep_lst) != len(head_lst)):
			print name,"Separator cannot match target completely!!"
			sys.exit(1)
		InOutput_Target=[]
		for i in range(0,len(tar_lst)):
			InOutput_Target.append(head_lst[i]+'=('+sep_lst[i].join(tar_lst[i])+')\n')
		InOutput_Target.append('\n')
		return InOutput_Target

class Input_Output_Settings(Read_ConfigFile,InOutput_Target_Pattern):
	def __init__(self, outpath4proj,module,step,sample,config):
        	self.outpath4proj = outpath4proj
        	self.sample = sample 
		self.step = step
		self.module=module
		Read_ConfigFile.__init__(self,config)
		InOutput_Target_Pattern.__init__(self,"")
	
	def Process_SampleInfo(self):
		if (re.findall('fq.gz|fastq.gz',self.sample)):
			(self.name,self.fq)=self.sample.split('|')
		elif(re.search('^FamGroup|^chr|all\||F..',self.sample)):
			(self.name,self.mem)=self.sample.split('|')
		else:
			self.name=self.sample
	
	def Define_Outdirs(self):
		self.outpath4step=self.outpath4proj+'/'+self.module+'/'+self.name+'/'+self.step
		if  not os.path.exists(self.outpath4step):
			os.makedirs(self.outpath4step)
	
	def Define_Softwares(self,sec,opt):
		self.software=self.Get_OptionLists4shell(sec,opt,0)

	def Define_Parameters(self,sec,opt):
		self.parameter=self.Get_OptionLists4shell(sec,opt,1)	
	
	def Define_References(self,sec,opt):
		self.reference=self.Get_OptionLists4shell(sec,opt,2)
	
	def Define_InHouse_Scripts(self,sec,opt):
		self.InHouse_Scripts=self.Get_OptionLists4shell(sec,opt,3)

	def Define_InputTargets(self,head_lst):
		self.InputTarget=self.Get_OptionItems(self.step + ".input","input.target")
		self.InputPath=self.Get_OptionItems(self.step + ".input","input.path")
		if (self.InputTarget=="from_sample.ini"):
			(R1,R2,self.InputPath)=self.fq.split(":")
			Mem=[R1,R2]
			Path_lst=['${InputPath}/','${InputPath}/']
			tar=self.Target_Simp(Path_lst,Mem,'\s+')
			sep_lst=[' ',' ',' ']
			tar_lst=[tar[0],tar[1],['${fq1[@]} ${fq2[@]}']]
			self.InputTarget=self.Target_Join(sep_lst,tar_lst,head_lst)
		elif(self.InputTarget=="from_fam.ped"):
			self.Pattern=self.Get_OptionItems(self.step + ".input","input.pattern")
			tar=self.Target_Lst(self.Pattern,self.mem,':')
			sep_lst=[' ']
			tar_lst=[tar[0]]
			self.InputTarget=self.Target_Join(sep_lst,tar_lst,head_lst)
		elif(self.InputTarget=="from_fam.group"):
			self.Pattern=self.Get_OptionItems(self.step + ".input","input.pattern")
			tar=self.Target_Lst(self.Pattern,self.mem,':')
			sep_lst=[' ',' -V ']
			tar_lst=[tar[0],tar[0]]
			self.InputTarget=self.Target_Join(sep_lst,tar_lst,head_lst)
		elif(self.InputTarget=="from_chr.group"):
			self.Pattern=self.Get_OptionItems(self.step + ".input","input.pattern")
			tar=self.Target_Lst(self.Pattern,self.mem,':')
			sep_lst=[' ',' -V ']
			tar_lst=[tar[0],tar[0]]
			self.InputTarget=self.Target_Join(sep_lst,tar_lst,head_lst)
		elif(self.InputTarget=="from_chrfam.group"):
			self.Pattern=self.Get_OptionItems(self.step + ".input","input.pattern")
			tar=self.Target_Lst(self.Pattern,self.mem,':')
			sep_lst=[' ',' --variant ']
			tar_lst=[tar[0],tar[0]]
			self.InputTarget=self.Target_Join(sep_lst,tar_lst,head_lst)
		elif(self.InputTarget == 'from_all.group'):
			self.Pattern=self.Get_OptionItems(self.step + ".input","input.pattern")
			tar=self.Target_Lst(self.Pattern,self.mem,':')
			sep_lst=['  ','  ','  ','  ','  ']
			tar_lst=[sum(tar,[]),tar[0],tar[1],tar[2],tar[3]]
			self.InputTarget=self.Target_Join(sep_lst,tar_lst,head_lst)
		else:
			Path_lst=['${InputPath}/']
			Mem_lst=[self.InputTarget]
			tar_lst=self.Target_Simp(Path_lst,Mem_lst,':')
			sep_lst=['  ']
			self.InputTarget=self.Target_Join(sep_lst,tar_lst,head_lst)

	def Define_OutputTargets(self,head_lst):
		self.OutputTarget=self.Get_OptionItems(self.step + ".output","output.target")
		self.OutputPath=self.Get_OptionItems(self.step + ".output","output.path")
		Path_lst=['${OutputPath}/']
		Mem_lst=[self.OutputTarget]
		tar_lst=self.Target_Simp(Path_lst,Mem_lst,':')
		sep_lst=[' ']
		self.OutputTarget=self.Target_Join(sep_lst,tar_lst,head_lst)
	 	
	def Define_Basic(self):
		self.basic=['###Basic\n','outpath4proj='+self.outpath4proj+'\n']
		self.basic.append('name='+self.name+'\n')
		self.basic.append('outpath4step='+self.outpath4step+'\n')
		self.basic.append('InputPath='+self.InputPath+'\n')
		self.basic.append('OutputPath='+self.OutputPath+'\n')
		self.basic.extend(self.InputTarget)
		self.basic.extend(self.OutputTarget)

def Define_Requirements(outpath4proj,module,step,sample,config,sec,opt):	
	requirements=Input_Output_Settings(outpath4proj,module,step,sample,config)
	requirements.Process_SampleInfo()
	requirements.Define_Outdirs()
	requirements.Define_Softwares(sec,opt)
	requirements.Define_Parameters(sec,opt)
	requirements.Define_References(sec,opt)
	requirements.Define_InHouse_Scripts(sec,opt)
	requirements.Define_InputTargets(opt[-2])
	requirements.Define_OutputTargets(opt[-1])
	requirements.Define_Basic()
	return requirements

def Check_InputOutputTarget(step):
	log1=['for file in ${inputtarget[@]}\n',
		'do\n',
        	'if [ ! -f "${file}" ]; then\n',
        	'\techo -e the target file ${file} cannot be generated successfully!!\n'
        	'\texit 0\n',
                'fi\n',
                'done\n',
		'echo -e The Step '+step+' started at `date` "\\n" >>${outpath4step}/run.log;\n\n']
	
        log2=[	'echo -e The Step '+step+' ended at `date` "\\n" >>${outpath4step}/run.log;\n',
		'for file in ${outtarget[@]}\n',
        	'do\n',
        	'if [ ! -f "${file}" ]; then\n',
        	'\techo -e the target file ${file} cannot bed generated successfully!!\n',
        	'\texit 0\n',
        	'fi\n',
        	'done\n',
        	'\tmv ${outpath4step}/run.log ${outpath4step}/'+step+'.log\n\n']
	return(log1,log2)

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

def Write_makefile(outfile,outpath,tar_out,ext_out,actStep):
	step_tar=['outpath=',outpath,'\n',
		'\nstep_tar=',tar_out,
		'\n.PHONY:all\n',
		'all:$(step_tar)\n']
	step0=['define steps\n']
	step1=['$2:'+'$(outpath)/$1/'+actStep+'/'+actStep+'.sh\n',
		'\tsh $(outpath)/$1/'+actStep+'/'+actStep+'.sh\n']
	step2=['endef\n']
	output=sum([step_tar,step0,step1,step2,ext_out],[])
	f=open(outfile+actStep+'.mk',"w")
	f.writelines(output)
	f.close

def Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands): 
        ### Basic info, software, parameter,reference
        sec=['software',step+'.parameter','reference','inhouse_scripts']
	requirements=Define_Requirements(outpath4proj,module,step,sample,config,sec,opt)
        require=requirements
        require_setting=[]
        require_setting=sum([require.basic,require.software,require.parameter,require.reference,require.InHouse_Scripts],[])
        (log1,log2)=Check_InputOutputTarget(step)
	### write makefile
	output=sum([require_setting,log1,commands,log2],[])
	makefilename=require.outpath4step+'/'+step+'.sh'
	f=open(makefilename,"w")
	output=sum([require_setting,log1,commands,log2],[])
	f.writelines(output)
	f.close

if __name__ == '__main__':
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
	
	#requirements preparation
	actSteps=[]
	ConfigInfo=Read_ConfigFile(config)
	module_path=ConfigInfo.Get_OptionItems('software','module_path')
	steps=ConfigInfo.Get_OptionItems('actSteps','StepstoRun')
	actSteps=steps.split(',')
	waitSteps=ConfigInfo.Get_OptionItems('actSteps','StepstoWait')
	print actSteps, waitSteps
