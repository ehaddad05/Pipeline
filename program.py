import pysam
import argparse
import os.path
import time
import pandas as pd
import sys
import argparse
import subprocess
import importlib
from functools import reduce
from SRAClass import SRA
from GBKClass import GBK
from ProjectClass import Project


from Bio import SeqIO



def parser():
	parser = argparse.ArgumentParser()
	### Args for processing an accession file
	parser.add_argument("-p",'--process', help="Option to process data", dest = "process")
	parser.add_argument("-a",'--accession', help="SRA accession file")
	parser.add_argument("-n",'--name', help="SRA accession file descriptive name")
	parser.add_argument("-g",'--gbk', help="GBK file")
	### Args for creating a project out of specific accession files
	parser.add_argument("-u",'--project', help = "Option to make a new project (provide name of project here too)",dest = "projectname")
	parser.add_argument("-q",'--quick', nargs='+', help = "Use this option if you are directly providing fcount references, and list the descriptive names of those references",dest="FCountnames")
	parser.add_argument("-l",'--list', nargs='+', help="<Required> Provide list of descriptive names of accessions or fcount references",dest="list")
	parser.add_argument("-k", '--kegg', help = "Provide kegg reference file here", dest = "kegg")
	## Args for only getting specific things
	parser.add_argument("-t", '--get', help = "Specifies what you want to get. Current Options: GBK amino acid fasta (gbk_fa), GBK nucleotide fasta (gbk_na)", dest = "get")
	parser.add_argument("-i", '--input', help = "Provide input here", dest = "input")
	args = parser.parse_args()


	return(args)


def CacheData(SRA,GBK):
	if not os.path.exists("/home/kz/Pipeline/cache.txt"):
		with open('cache.txt', 'w') as f:
			organization = { "SRA" : {}, "GBK" : {}}
			SRAObject = SRA.OrganizeData()
			GBKObject = GBK.OrganizeData()
			organization["SRA"].update({SRA.DescriptiveName : {} })
			organization["GBK"].update({GBK.LocusTag : {} })
			organization["SRA"][SRA.DescriptiveName].update(SRAObject)
			organization["GBK"][GBK.LocusTag].update(GBKObject)
			print(organization,file=f)
	else:
		with open('cache.txt', 'r+') as f:
			content = f.read();
			organization = eval(content)
			SRAObject = SRA.OrganizeData()
			GBKObject = GBK.OrganizeData()
			organization["SRA"].update({SRA.DescriptiveName : {} })
			organization["GBK"].update({GBK.LocusTag : {} })
			organization["SRA"][SRA.DescriptiveName].update(SRAObject)
			organization["GBK"][GBK.LocusTag].update(GBKObject)
			f.truncate(0)
			print(organization,file=f)		
    		

def main(args):

	if hasattr(args,'process') and args.process is not None:

		InputSRA = SRA(args.n, args.a)
		print("SRA is downloaded")
		InputGBK = GBK(args.g)
		print("GBK is read")
		InputGBK.MakeSAF()
		InputGBK.MakeNa("MIDDLE")
		InputGBK.MakeFa("MIDDLE")
		print("GBK SAF is made")
		InputSRA.TrimAccession()
		print("SRA is trimmed")
		InputSRA.BuildIndex(InputGBK)
		print("Index is built")
		InputSRA.bowtie2_align(InputGBK)
		print("SRA Outputs retrieved")
		InputSRA.FeatureCounts(InputGBK)
		print("Featurecounts finished")
		CacheData(InputSRA,InputGBK)

	elif hasattr(args,'projectname') and args.projectname is not None:
		if not hasattr(args, 'FCountnames'):
			with open('cache.txt', 'r+') as f:
				content = f.read();
				organization = eval(content)
			SRATable = {}
			for i in range(0,len(args.list)):
				SRATable[i] = (SRA(organization["SRA"][args.list[i]]))
			NewProject = Project(args.projectname,SRATable, "SRA")
		else:
			fcounts = {}
			for i in range(0,len(args.list)):
				fcounts[args.FCountnames[i]] = args.list[i]
			NewProject = Project(args.projectname,fcounts, "Fcount",args.kegg)
		NewProject.MakeFCountDataFrame()
		NewProject.NormalizeReads()
		NewProject.CalculateTE()
		NewProject.MergeWithKegg()
	elif hasattr(args,'get') and args.get is not None:
		InputGBK = GBK(args.input)
		if args.get == "gbk_fa":
			InputGBK.MakeFa("OUTPUT")
		elif args.get == "gbk_na":
			InputGBK.MakeNa("OUTPUT")
		else:
			print("No \"get\" option was selected")
		
	else:
		print("No option was selected")


	


if __name__ == "__main__":
    args = parser()
    main(args)




















