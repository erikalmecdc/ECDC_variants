#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
from glob import glob
from datetime import datetime
import sys

#The current variant table (refer to https://www.ecdc.europa.eu/en/covid-19/variants-concern for timestamp)
url='https://www.ecdc.europa.eu/sites/default/files/documents/PathogenVariant_public_mappings.csv'

url_static = 'https://raw.githubusercontent.com/erikalmecdc/ECDC_variants/main/data/PathogenVariant_info.csv'
url_static_VUM = 'https://raw.githubusercontent.com/erikalmecdc/ECDC_variants/main/data/PathogenVariant_info_VUM.csv'


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

def VariantsQuery(df_ECDC_variants_in,query_pango_lineage=None,query_mutations=None,df_ECDC_variants_VUM_in=pd.DataFrame()):

	ECDCVariant = 'Not classified by ECDC'
	ECDCVariants_VUM = []

	for VirusVariant in df_ECDC_variants_in.index:
		
		Sublineages = df_ECDC_variants_in.loc[VirusVariant,'included sub-lineages'].split('|')
		ECDCClassification = df_ECDC_variants_in.loc[VirusVariant,'ECDCClassification']
		LineageMutations = {}

		if not pd.isnull(df_ECDC_variants_in.loc[VirusVariant,'LineageMutations']):
			LineageMutations =  json.loads(df_ECDC_variants_in.loc[VirusVariant,'LineageMutations'])


		if query_pango_lineage:

			if query_pango_lineage in Sublineages:
				ECDCVariant = VirusVariant+' ('+ECDCClassification+')'

		if query_mutations and LineageMutations:

			query_mutations = set([x.strip(' ') for x in query_mutations.split(',')])

			if not query_pango_lineage:
				query_pango_lineage = 'any_pango_lineage'

			if query_pango_lineage in LineageMutations.keys():
				required_mutations = set(LineageMutations[query_pango_lineage].split('+'))
				
				if len(required_mutations - query_mutations) == 0:
					ECDCVariant = VirusVariant+' ('+ECDCClassification+')'

	if not df_ECDC_variants_VUM_in.empty:

		for VirusVariant in df_ECDC_variants_in.index:
		
			Sublineages = df_ECDC_variants_in.loc[VirusVariant,'included sub-lineages'].split('|')
			ECDCClassification = df_ECDC_variants_in.loc[VirusVariant,'ECDCClassification']
			LineageMutations = {}

			if not pd.isnull(df_ECDC_variants_in.loc[VirusVariant,'LineageMutations']):
				LineageMutations =  json.loads(df_ECDC_variants_in.loc[VirusVariant,'LineageMutations'])


			if query_pango_lineage:

				if query_pango_lineage in Sublineages:
					ECDCVariant_VUM.append(VirusVariant+' ('+ECDCClassification+')')

			if query_mutations and LineageMutations:

				query_mutations = set([x.strip(' ') for x in query_mutations.split(',')])

				if not query_pango_lineage:
					query_pango_lineage = 'any_pango_lineage'

				if query_pango_lineage in LineageMutations.keys():
					required_mutations = set(LineageMutations[query_pango_lineage].split('+'))
				
					if len(required_mutations - query_mutations) == 0:
						ECDCVariant_VUM.append(VirusVariant+' ('+ECDCClassification+')')

	if len(ECDCVariants_VUM) == 0:
		ECDCVariants_VUM = 'Not classified by ECDC'
	else:
		'/'.join(ECDCVariants_VUM)

	return ECDCVariant,ECDCVariants_VUM

def EpiCoVmetadataQuery(df_in,df_ECDC_variants_in):

	df_in['ECDCClassification'] = ''
	df_in['ECDCClassification_VUM'] = ''

	for index in df_in.index:

		Pango_lineage = df_in.loc[index,'Lineage'].split(' ')[0]
		AA_Substitutions = df_in.loc[index,'AA Substitutions'].replace('(','').replace(')','')
		ECDCClassification , ECDCClassification_VUM = VariantsQuery(df_ECDC_variants_in,Pango_lineage,AA_Substitutions)

		df_in.loc[index,'ECDCClassification'] = ECDCClassification
		df_in.loc[index,'ECDCClassification_VUM'] = ECDCClassification_VUM

	return df_in




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


	df_ECDC_variants=pd.read_csv(url_static,index_col='VirusVariant')
	df_ECDC_variants_VUM=pd.read_csv(url_static,index_col='VirusVariant')

	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants,'Q.7')
	print(ECDCVariant+'\n')

	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants,'HK.22')
	print(ECDCVariant+'\n')

	#pango lineage with comma separated mutations that does not meet requirement
	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants,'FD.1','Spike_F486P, Spike_F490S',df_ECDC_variants_VUM)
	print(ECDCVariant+'\n')

	#pango lineage with comma separated mutations that does not meet requirement
	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants,'XBB.1.5','Spike_L455F, Spike_F456L',df_ECDC_variants_VUM)
	print(ECDCVariant+'\n')
	print(ECDCVariantsVUM+'\n')

	#pango lineage with comma separated mutations that meets requirement
	ECDCVariant,ECDCVariantsVUM = VariantsQuery(df_ECDC_variants,'FD.1','Spike_F486P, Spike_F456L, Spike_F490S')
	print(ECDCVariant+'\n')

	df_EpiCoVmetadata,snapshot_date = readInputFromEpiCoVmetadataFile()

	if not df_EpiCoVmetadata.empty:
		df_EpiCoVmetadata = EpiCoVmetadataQuery(df_EpiCoVmetadata,df_ECDC_variants)
		outfile_EpiCoVmetadata = 'output/EpiCoVmetadata_'+str(snapshot_date)+'.xlsx'
		df_EpiCoVmetadata.to_excel(outfile_EpiCoVmetadata)
		print('ECDC variants assigned to '+outfile_EpiCoVmetadata)