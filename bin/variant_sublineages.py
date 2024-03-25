#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
from glob import glob
from datetime import datetime
import sys

#The current variant table (refer to https://www.ecdc.europa.eu/en/covid-19/variants-concern for timestamp)
url='https://www.ecdc.europa.eu/sites/default/files/documents/PathogenVariant_public_mappings.csv'
url_VUM = 'https://www.ecdc.europa.eu/sites/default/files/documents/PathogenVariant_public_mappings_VUM.csv'

#These are similar files on git for testing purposes
url_static = 'https://raw.githubusercontent.com/erikalmecdc/ECDC_variants/main/data/PathogenVariant_info.csv'
url_static_VUM = 'https://raw.githubusercontent.com/erikalmecdc/ECDC_variants/main/data/PathogenVariant_public_mappings_VUM.csv'

EpiPulse_template = pd.read_csv('data/NCOVVARIANT_template.csv',index_col='RecordId')

limit = 3

#headers of table:
#VirusVariant => The variant's name referred to by ECDC.
#Sublineages => All Pangolin lineages that are included in variant independent of amino acid substitutions.
#LineageMutations => json text that, converted to a dictionary, has Pangolin lineage as key and required amino acid substitutions as values. The amino acid substitutions are separated by +. 
#ECDCClassification => Classification of the variant, i.e. VOC (variant of concern), VOI (variant of interest) or VUM (variant under monitoring).

#Describes the content of the table in detail from examples
def DescribeVariantTable(df_in):
	for VirusVariant in df_in.index:

		Sublineages = df_in.loc[VirusVariant,'included sub-lineages'].split('|')
		ECDCClassification = df_in.loc[VirusVariant,'ECDCClassification']
		
		LineageMutations = {}

		if not pd.isnull(df_in.loc[VirusVariant,'LineageMutations']):
			LineageMutations =  json.loads(df_in.loc[VirusVariant,'LineageMutations'])

		LineageMutationsExamples = ''
		if LineageMutations:
			i = 0

			LineageMutationsExamples += ' for '+str(len(LineageMutations.keys()))+' sub-lineages e.g.'
			for lineage, mutations in LineageMutations.items():
				mutations = mutations.split('+')

				lineage = lineage.replace('any_pango_lineage','any pango lineage')
				LineageMutationsExamples+= ' '+lineage+' with '+'+'.join(mutations)+','

				i += 1

				if i >= limit:
					LineageMutationsExamples = LineageMutationsExamples.rstrip(',')
					break

		VariantInfo = '\nVariant '+VirusVariant+ ' is classified by ECDC as a '+ECDCClassification+' and includes '+str(len(Sublineages))+ ' sub-lineages such as '+','.join(Sublineages[0:limit])

		if LineageMutationsExamples:

			VariantInfo += ',\nwith mutation criterium'+LineageMutationsExamples

		VariantInfo += '.'

		print(VariantInfo)

#Assigns ECDC status for VOI/VOC (also selected de-escalated) and VUM(s)
def VariantsQuery(df_ECDC_variants_in,query_pango_lineage=None,query_mutations=None,df_ECDC_variants_VUM_in=pd.DataFrame(),recordid=None):

	ECDCVariant = 'Not classified by ECDC'
	ECDCVariants_VUM = []

	VirusVariantDetected = ''

	for VirusVariant in df_ECDC_variants_in.index:
		
		Sublineages = df_ECDC_variants_in.loc[VirusVariant,'included sub-lineages'].split('|')
		ECDCClassification = df_ECDC_variants_in.loc[VirusVariant,'ECDCClassification']
		LineageMutations = {}

		if not pd.isnull(df_ECDC_variants_in.loc[VirusVariant,'LineageMutations']):
			LineageMutations =  json.loads(df_ECDC_variants_in.loc[VirusVariant,'LineageMutations'])


		if query_pango_lineage:

			if query_pango_lineage in Sublineages:
				ECDCVariant = VirusVariant+' ('+ECDCClassification+')'
				VirusVariantDetected = VirusVariant

		if query_mutations and LineageMutations:

			query_mutations = set([x.strip(' ') for x in query_mutations.split(',')])

			if not query_pango_lineage:
				query_pango_lineage = 'any_pango_lineage'

			if query_pango_lineage in LineageMutations.keys():
				required_mutations = set(LineageMutations[query_pango_lineage].split('+'))
				
				if len(required_mutations - query_mutations) == 0:
					ECDCVariant = VirusVariant+' ('+ECDCClassification+')'
					VirusVariantDetected = VirusVariant

	if recordid and ECDCVariant != 'Not classified by ECDC':
		EpiPulse_template.loc[recordid,'VirusVariant'] = VirusVariantDetected
	elif recordid and query_pango_lineage:
		EpiPulse_template.loc[recordid,'VirusVariantOther'] = query_pango_lineage

	if not df_ECDC_variants_VUM_in.empty:

		for VirusVariant in df_ECDC_variants_VUM_in.index:
		
			Sublineages = df_ECDC_variants_VUM_in.loc[VirusVariant,'included sub-lineages'].split('|')
			ECDCClassification = df_ECDC_variants_VUM_in.loc[VirusVariant,'ECDCClassification']
			LineageMutations = {}

			if not pd.isnull(df_ECDC_variants_VUM_in.loc[VirusVariant,'LineageMutations']):
				LineageMutations =  json.loads(df_ECDC_variants_VUM_in.loc[VirusVariant,'LineageMutations'])


			if query_pango_lineage:

				if query_pango_lineage in Sublineages:
					ECDCVariants_VUM.append(VirusVariant+' ('+ECDCClassification+')')

			if query_mutations and LineageMutations:

				query_mutations = set([x.strip(' ') for x in query_mutations.split(',')])

				if not query_pango_lineage:
					query_pango_lineage = 'any_pango_lineage'

				if query_pango_lineage in LineageMutations.keys():
					required_mutations = set(LineageMutations[query_pango_lineage].split('+'))
				
					if len(required_mutations - query_mutations) == 0:
						ECDCVariants_VUM.append(VirusVariant+' ('+ECDCClassification+')')
	else:
		ECDCVariants_VUM = ['VUM status not examined']

	if len(ECDCVariants_VUM) == 0:
		ECDCVariants_VUM = 'Not classified as a VUM by ECDC'
	else:
		ECDCVariants_VUM = '/'.join(ECDCVariants_VUM)

	return ECDCVariant,ECDCVariants_VUM

#Multi assignment of ECDC classifications by GISAID EpiCoV patient status metadata input, file should have this format gisaid_hcov-19_2024_03_22_08.tsv
def EpiCoVmetadataQuery(df_in,df_ECDC_variants_in):

	df_in['ECDCClassification'] = ''
	df_in['ECDCClassification_VUM'] = ''

	for index in df_in.index:

		Pango_lineage = df_in.loc[index,'Lineage'].split(' ')[0]
		AA_Substitutions = df_in.loc[index,'AA Substitutions'].replace('(','').replace(')','')
		ECDCClassification , ECDCClassification_VUM = VariantsQuery(df_ECDC_variants_in,Pango_lineage,AA_Substitutions,recordid=index)

		df_in.loc[index,'ECDCClassification'] = ECDCClassification
		df_in.loc[index,'ECDCClassification_VUM'] = ECDCClassification_VUM

	return df_in

#Reading GISAID EpiCoV patient metadata from input folder
def readInputFromEpiCoVmetadataFile():

	EpiCoVmetadataFiles = sorted(glob('input/*.tsv'))

	if len(EpiCoVmetadataFiles) > 0:
		
		print('Found EpiCoV metadata file: '+EpiCoVmetadataFiles[-1]+'\n')
		
		snapshot_date_out = '-'.join(EpiCoVmetadataFiles[-1].split('_')[2:5]).replace('.tsv','')
		snapshot_date_out = datetime.strptime(snapshot_date_out, "%Y-%m-%d").date()
		#print(snapshot_date_out)

		return  pd.read_csv(EpiCoVmetadataFiles[-1], sep='\t',index_col='Virus name'),snapshot_date_out

	else:
		return pd.DataFrame(),''


if __name__ == "__main__":

	#Collect variant statuses (Switch to the non-static variables to use up to date data)
	df_ECDC_variants=pd.read_csv(url_static,index_col='VirusVariant')
	df_ECDC_variants_VUM=pd.read_csv(url_static_VUM,index_col='VirusVariant')

	#Run this to get a detailed description of the tables
	DescribeVariantTable(df_ECDC_variants)

	#A single query of a de-escalated VOC, both VOI/VOC and VUM status
	#Also a recordid is provided and activate addition to an EpiPulse/TESSy reporting sheet
	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants_in=df_ECDC_variants,query_pango_lineage='Q.7',df_ECDC_variants_VUM_in=df_ECDC_variants_VUM,recordid='Domestic_sample1')
	print(ECDCVariant)
	print(ECDCVariantsVUM+'\n')

	#A single query for a variant monitored as a VOI during the production of the test data. Since VUM table is not provided this status is not examined
	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants_in=df_ECDC_variants,query_pango_lineage='HK.22')
	print(ECDCVariant)
	print(ECDCVariantsVUM+'\n')

	#A single query for using mutations that are not defining any variant
	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants_in=df_ECDC_variants,query_pango_lineage='FD.1',query_mutations='Spike_F486P, Spike_F490S')
	print(ECDCVariant+'\n')

	#A single query for using mutations that are defining the VOI XBB.1.5-like+F456L
	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants_in=df_ECDC_variants,query_pango_lineage='XBB.1.5',query_mutations='Spike_L455F, Spike_F456L')
	print(ECDCVariant+'\n')

	#A single query for a currently monitored VUM, notice that it counts as a de-escalated VOI BA.2 parent as well
	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants_in=df_ECDC_variants,query_pango_lineage='BA.2.87.1',df_ECDC_variants_VUM_in=df_ECDC_variants_VUM)
	print(ECDCVariant)
	print(ECDCVariantsVUM+'\n')


	#Instead of processing one entry, process arbitrary GISAID EpiCov metadata
	df_EpiCoVmetadata,snapshot_date = readInputFromEpiCoVmetadataFile()

	if not df_EpiCoVmetadata.empty:
		df_EpiCoVmetadata = EpiCoVmetadataQuery(df_EpiCoVmetadata,df_ECDC_variants)
		outfile_EpiCoVmetadata = 'output/EpiCoVmetadata_'+str(snapshot_date)+'.xlsx'
		df_EpiCoVmetadata.to_excel(outfile_EpiCoVmetadata)
		print('ECDC variants assigned to '+outfile_EpiCoVmetadata)

	if not EpiPulse_template.empty:
		outfile_EpiPulse = 'output/EpiPulse_'+str(snapshot_date)+'.xlsx'
		EpiPulse_template.to_excel(outfile_EpiPulse)
		print('An EpiPulse/TESSy import has been written to '+outfile_EpiPulse)