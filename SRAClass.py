import pysam
import argparse
import os.path
import time
import pandas as pd
import argparse
import sys
import subprocess
from Bio import SeqIO		

def cmd(cmd, path):  
    return subprocess.check_output(cmd, cwd=path, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')

def Command(*args, **kwargs):
    fullString = ""
    currentDirectory = ""
    for a in args:
        fullString+= (str(a) + " ")
    for k,v in kwargs.items():
        currentDirectory = v
    print("cmd: %s" % fullString)
    cmd(fullString,currentDirectory);  


def GenerateLink(Accession):
	###ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR769/005/SRR7698518/SRR7698518.fastq.gz
	LinkStart = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
	LinkEnd = ".fastq.gz"
	LinkGen = Accession[0:6] + "/00" + Accession[-1:] + "/" +Accession + "/" + Accession
	return(LinkStart+LinkGen+LinkEnd)
	




class SRA:

	def __init__(self,DescriptiveName, Accession):
		self.DescriptiveName = DescriptiveName
		self.Accession = Accession
		if not os.path.exists("/home/kz/Pipeline/SRA_INPUT_FILES/%s.fastq" % self.Accession) and not os.path.exists("/home/kz/Pipeline/SRA_INPUT_FILES/%s.fastq.gz" % self.Accession):
			try:
				Command("wget ", GenerateLink(Accession), Path = "/home/kz/Pipeline/SRA_INPUT_FILES")
				self.Reference = "/home/kz/Pipeline/SRA_INPUT_FILES/" + Accession + ".fastq.gz"
			except:
				Command("fastq-dump ", self.Accession, Path = "/home/kz/Pipeline/SRA_INPUT_FILES")
				self.Reference = "/home/kz/Pipeline/SRA_INPUT_FILES/" + Accession + ".fastq"
		
			while not os.path.exists(self.Reference):
				time.sleep(1)
			if self.Reference[-3] == ".gz":
				Command("gunzip", self.Reference, Path = "/home/kz/Pipeline/SRA_INPUT_FILES")
				self.Reference = "/home/kz/Pipeline/SRA_INPUT_FILES/" + Accession + ".fastq"
				while not os.path.exists(self.Reference):
					time.sleep(1)
		self.Reference = "/home/kz/Pipeline/SRA_INPUT_FILES/" + Accession + ".fastq"

	def __init__(self,cacheDict):
		self.DescriptiveName = cacheDict["DescriptiveName"]
		self.Accession = cacheDict["Accession"]
		self.Reference = cacheDict["Reference"] 
		self.TrimmedReference = cacheDict["TrimmedReference"] 
		self.SAMReference = cacheDict["SAMReference"] 
		self.BOWReference = cacheDict["BOWReference"]
		self.BAMReference = cacheDict["BAMReference"]
		self.sortedBAMReference = cacheDict["sortedBAMReference"] 
		self.FcountReference = cacheDict["FcountReference"] 	

	def OrganizeData(self):
		SRAObject = {"DescriptiveName" : self.DescriptiveName, "Accession" : self.Accession, "Reference" : self.Reference, "TrimmedReference" : self.TrimmedReference, "SAMReference" : self.SAMReference, "BOWReference" : self.BOWReference, "BAMReference" : self.BAMReference, "sortedBAMReference" : self.sortedBAMReference, "FcountReference" : self.FcountReference}
		return SRAObject

	def TrimAccession(self):
		self.TrimmedReference = "/home/kz/Pipeline/TRIMMED_SRA_INPUT_FILES/" + self.Accession +"_trimmed.fq"
		if not os.path.exists(self.TrimmedReference):
			Command("trim_galore -fastqc ", self.Reference, Path = "/home/kz/Pipeline/TRIMMED_SRA_INPUT_FILES")
		
		
	
	def BuildIndex(self, GBK):
		Command("bowtie2-build ", GBK.FastaReference, " " + GBK.IndexName, Pipe = "/home/kz/Pipeline/MIDDLE/")
		return self	
	def bowtie2_align(self,GBK):
		self.SAMReference = "/home/kz/Pipeline/MIDDLE/" + self.Accession + ".sam"
		self.BOWReference = "/home/kz/Pipeline/MIDDLE/" + self.Accession + ".bow"
		Command("bowtie2", " -x ", GBK.IndexName, " -U ", self.TrimmedReference, " -p 8 ", " -S ", self.SAMReference, " 2> ", self.Accession + ".bow", Pipe = "/home/kz/Pipeline/MIDDLE/")
		self.sortedBAMReference = "/home/kz/Pipeline/MIDDLE/sorted_" + self.Accession + "_1.bam"
		self.BAMReference = "/home/kz/Pipeline/MIDDLE/" + self.Accession + ".bam"
		Command("samtools view -S -b ", self.SAMReference, " > ", self.Accession + ".bam", Pipe = "/home/kz/Pipeline/MIDDLE/")
		bamfile = pysam.AlignmentFile(self.BAMReference, "rb") 
		pysam.sort("-o", self.sortedBAMReference, self.BAMReference)

		
	

	def FeatureCounts(self,GBK):
		self.FcountReference = "/home/kz/Pipeline/OUTPUT/%s.fcount" % self.Accession
		Command("featureCounts", "-a %s" % GBK.SAF, " -F SAF", "-o %s.fcount" % self.Accession, "-P -C -d 30 -D 800", "%s" % self.sortedBAMReference, Pipe = "/home/kz/Pipeline/OUTPUT/")
	
				
	def ChangeName(self,Name):
		self.DescriptiveName = Name

