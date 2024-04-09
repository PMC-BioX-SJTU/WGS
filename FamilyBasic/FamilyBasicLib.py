#/usr/bin/python
import ConfigParser
import string
import os
import argparse
import sys
import re
sys.path.append('/home/zhujh/Codes/Pipelines/BasicWESorWGSPlus/')
from BasicWESorWGSLib import *

class Obtain_FamilyInfo():
        def __init__(self,pedFile):
                F=open(pedFile,'r')
		self.ped=F.readlines()
		F.close()
	
	def Get_famInfo(self,index):
		fam={}
		for i in range (0,len(self.ped)):
			line=self.ped[i].split()
			if(index == ':'):
				lines=[self.ped[i]]
			else:
				lines=[line[1]]
			if(line[0] not in fam):
				fam[ line[0] ]=lines
			else:
				fam[line[0]]+=lines
		return fam

	def Get_famMem(self):
		self.famMem=self.Get_famInfo(1)

	def Get_famGroup(self,count):
		self.famGroup=[]
		j=count
		lst=[]
		for i in range (0,len(self.ped)):
			line=self.ped[i].split()
			if(j>0):
				lst.append(line[1])
				j=j-1
			else:
				self.famGroup.append(lst)
				j=count
				lst=[line[1]]
				j=j-1
		self.famGroup.append(lst)
	
	def Output_famPed(self,outpath4step):
		famPed=self.Get_famInfo(':')
		for fam in famPed.iterkeys():
			outfile=outpath4step+'/'+fam+'/'+fam+'.ped'
			if  not os.path.exists(outpath4step+'/'+fam):
				os.makedirs(outpath4step+'/'+fam)
			O=open(outfile,'w');
			O.writelines(famPed[fam])
			O.close()

	def Output_famPed_Trio(self,outpath4step):
		famPed=self.Get_famInfo(':')
		for fam in famPed.iterkeys():
			if ( len(famPed[fam]) == 3):
				outfile=outpath4step+'/'+fam+'.ped'
				O=open(outfile,'w');
				O.writelines(famPed[fam])
				O.close()
			else:
				k=0
				out=[]
				for i in range(0,len(famPed[fam])):
					items=famPed[fam][i].split()
					if items[2]!='0':
						k=k+1
						outfile=outpath4step+'/'+fam+'ch'+str(k)+'.ped'
						out.append(famPed[fam][i])
						O=open(outfile,'w')
						O.writelines(out)
						O.close()
						out.pop()
					else:
						out.append(famPed[fam][i])
						
	

if __name__ == "__main__" :
	ss=Obtain_FamilyInfo('../all.ped','BreakDancerCalling')
	ss.Get_famMem()
	ss.Define_InputTargets()
