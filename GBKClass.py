
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




class GBK:
	def __init__(self, Path):
		self.Reference = Path
		
		with open(self.Reference) as f:
			first_line = f.readline()
		Pieces = str.split(first_line, " ")
		
		self.LocusTag = Pieces[7]
		self.GenomeSize = Pieces[14]
		self.IndexName = self.LocusTag + "_index"
		
		



	def OrganizeData(self):
		GBKObject = {"Reference" : self.Reference, "LocusTag" : self.LocusTag, "GenomeSize" : self.GenomeSize, "IndexName" : self.IndexName, "ProteinReference" : self.ProteinReference, "FastaReference" : self.FastaReference, "SAF" : self.SAF}
		return GBKObject


	def MakeNa(self,where):
		self.FastaReference = "/home/kz/Pipeline/%s/GBK_na.fasta" % where
		with open(self.Reference, "r") as input_handle:
			with open(self.FastaReference, "w") as output_handle:
				sequences = SeqIO.parse(input_handle, "genbank")
				count = SeqIO.write(sequences, output_handle, "fasta")
	def MakeFa(self,where):
		cds=[]
		locus_tag=[]
		self.ProteinReference = "/home/kz/Pipeline/%s/%s.fa" % (where, self.LocusTag)

		rec = next(SeqIO.parse(self.Reference,'gb'))
		seq = rec.seq

		for feat in rec.features:
			if feat.type =='CDS':
				try:
					cds.append(feat.qualifiers['translation'][0])
				except:
					cds.append('na')
				try:   
					locus_tag.append(feat.qualifiers['locus_tag'][0])
				except:
					locus_tag.append('na')

		with open(self.ProteinReference,'w') as f:
			for i,j in zip(locus_tag,cds):
				f.write('>%s\n%s\n'%(i,j))

	def MakeSAF(self):
		self.SAF = "/home/kz/Pipeline/MIDDLE/%s.saf" % self.LocusTag
		infile = next(SeqIO.parse(self.Reference,'gb'))
		
		
		genes =[]
		strand = []
		start = []
		stop = []
		accession=infile.annotations['accessions'][0]

		for feature in infile.features:
		    if feature.type == 'CDS':  #Only protein-coding genes
		        genes.append(feature.qualifiers['locus_tag'][0])
		        if feature.strand == 1:
		            strand.append("+")
		            start.append(feature.location.start.real+1)
		            stop.append(feature.location.end.real)
		        elif feature.strand == -1:
		            strand.append("-")
		            start.append(feature.location.start.real+1)
		            stop.append(feature.location.end.real)
		gene_df = pd.DataFrame({"GeneID": genes,  'Chr':accession ,"Start": start, "Stop": stop,'Strand': strand}, columns = ['GeneID', 'Chr', 'Start','Stop', 'Strand'])
		gene_df.to_csv(self.SAF, sep='\t', index=False)

