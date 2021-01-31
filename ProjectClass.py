import pysam
import argparse
import os.path
import time
import pandas as pd
import json,sys,argparse
import subprocess
from functools import reduce
from Bio import SeqIO	



class Project:
	def __init__(self, Name, Table,t,kegg):
		self.Name = Name
		ReferenceList = {}
		DescriptiveList = {}
		if t == "SRA":
			Tbl = list(Table.values())
			for i in range(0,len(Tbl)):
				ReferenceList[i] = (Tbl[i].FcountReference)
				DescriptiveList[i] = (Tbl[i].DescriptiveName)
			self.ReferenceList = list(ReferenceList.values())
			self.DescriptiveList = list(DescriptiveList.values())
		elif t == "Fcount":
			ReferenceList = list(Table.values())
			DescriptiveList = list(Table.keys())
			self.ReferenceList = ReferenceList
			self.DescriptiveList = DescriptiveList
		self.KEGG = kegg
		Protein = []
		KEGG = []
		with open(self.KEGG) as myfile:
			for line in myfile:
				info = line.split('\t')
				Protein.append(info[0].strip())
				try:
					KEGG.append(info[1].strip())
				except:
					KEGG.append('na')
		self.KeggTable = "/home/kz/Pipeline/OUTPUT/"+self.Name+"_KEGG.csv"
		kegg_reads = pd.DataFrame({'Protein' : Protein, 'kegg' : KEGG})
		kegg_classification = pd.read_csv('/home/kz/Pipeline/kegg_classification_grouped_by_CAT_PATH.csv')
		kegg_merged = pd.merge(kegg_reads,kegg_classification,on=['kegg'])
		kegg_merged.to_csv(self.KeggTable, sep=",")
		
	
	def MakeFCountDataFrame(self):
		FCOUNT_DATA = []
		FCOUNT_INPUT = self.ReferenceList
		

		for FCOUNT in FCOUNT_INPUT:
			Protein = []
			Size = []
			Reads = []
			with open(FCOUNT) as myfile: 
				for line in myfile:
					if (len(line.split('\t')) == 1):
						continue
					info = line.split('\t')
					Protein.append(info[0])
					print(info[0])
					try:
						Size.append(int(info[5]))
					except:
						Size.append(info[5])
					try:
						Reads.append(int(info[6].strip()))
					except:
						Reads.append(info[6].strip())
			del Protein[0]
			del Size[0]
			del Reads[0]
			name = 'Reads_'+FCOUNT
			if len(self.DescriptiveList) > 0:
				name = self.DescriptiveList[FCOUNT_INPUT.index(FCOUNT)]
			FCOUNT_DATA.append(pd.DataFrame({'Protein' : Protein, 'gene_length' : Size, name : Reads}))
		self.TabulatedReference = "/home/kz/Pipeline/OUTPUT/"+self.Name+".csv"
		df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['Protein','gene_length'],how='inner'), FCOUNT_DATA)
		df_merged.to_csv(self.TabulatedReference, sep=",")

	def NormalizeReads(self):
		Reads = pd.read_csv(self.TabulatedReference)
		Reads = Reads.iloc[0:,1:]
		Columns = list(Reads.columns.values)
		Totals = {}
		ToDrop = []
		Found = False
		for i in range(0,len(Columns)):
			if Columns[i].find("Part") != -1:
				Found = True
				ToDrop.append(Columns[i])
				Piece = Columns[i][0:-5]
				if not Piece in Totals.keys():
					Totals[Piece] = []
				Totals[Piece].append(Columns[i])
		if Found:
			for key, tbl in Totals:
				for piece in tbl:
					Reads[key] += Reads[piece]
		
			CombinedReads = Reads.drop(columns=[ToDrop])
			CombinedReadsReference = "/home/kz/Pipeline/OUTPUT/combined_"+self.Name+".csv"
			CombinedReads.to_csv(self.CombinedReadsReference, sep=",")
			Reads = CombinedReads
		for i in range(0,len(Columns)):
			if Columns[i] != "Protein" and Columns[i] != "gene_length":
				Reads[Columns[i]+"_rpkm"]= (Reads[Columns[i]] * 1000000000)  /  (sum(Reads[Columns[i]]) * Reads["gene_length"])
			
		self.NormalizedReference = "/home/kz/Pipeline/OUTPUT/normalized_"+self.Name+".csv"
		Reads.to_csv(self.NormalizedReference, sep=",")
	def CalculateTE(self):
		Reads = pd.read_csv(self.NormalizedReference)
		Reads = Reads.iloc[0:,1:]
		Columns = list(Reads.columns.values)
		ProjectConditions = []
		Pairs = {}
		for i in range(0,len(Columns)):
			if Columns[i].find("rpkm") != -1:
				Condition = Columns[i][4:]
				ProjectConditions.append(Condition)
				
				if not Condition in Pairs.keys():
					Pairs[Condition] = {"RPF" : "s", "RNA" : "s"}
				if Columns[i].find("RPF") != -1: 
					
					Pairs[Condition]["RPF"] = Columns[i]
				
				elif Columns[i].find("RNA") != -1:
					
					Pairs[Condition]["RNA"] = Columns[i]
				
		for i in range(0,len(ProjectConditions)):
			c = ProjectConditions[i]
			Reads["TE_"+c[:-5]] = Reads[Pairs[c]["RPF"]] / Reads[Pairs[c]["RNA"]] 
		self.TECalcReference = "/home/kz/Pipeline/OUTPUT/TE_"+self.Name+".csv"
		Reads.to_csv(self.TECalcReference, sep=",")
	def MergeWithKegg(self):
		Reads = pd.read_csv(self.TECalcReference)
		Reads = Reads.iloc[0:,1:]
		Kegg = pd.read_csv(self.KeggTable)
		Kegg = Kegg.iloc[0:,1:]
		Final = pd.merge(Kegg, Reads,how="outer", on=['Protein'])
		self.MergedWithKegg = "/home/kz/Pipeline/OUTPUT/TE_wKegg_"+self.Name+".csv"
		Final.to_csv(self.MergedWithKegg, sep=",")



