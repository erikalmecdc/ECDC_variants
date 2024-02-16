#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json

#The current variant table (refer to https://www.ecdc.europa.eu/en/covid-19/variants-concern for timestamp)
url='https://www.ecdc.europa.eu/sites/default/files/documents/PathogenVariant_public_mappings.csv'

url_static = 'https://raw.githubusercontent.com/erikalmecdc/ECDC_variants/main/data/PathogenVariant_info.csv'

#An example of a variant table that is not maintained
filepath = 'df_PathogenVariant_info_CO_public.csv'

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

def SingleQuery(df_in,query_pango_lineage=None,query_mutations=None,VOCVOI=True):

	ECDCVariants = []

	for VirusVariant in df_in.index:
		
		Sublineages = df_in.loc[VirusVariant,'included sub-lineages'].split('|')
		ECDCClassification = df_in.loc[VirusVariant,'ECDCClassification']
		LineageMutations = {}

		if not pd.isnull(df_in.loc[VirusVariant,'LineageMutations']):
			LineageMutations =  json.loads(df_in.loc[VirusVariant,'LineageMutations'])


		if query_pango_lineage:

			if query_pango_lineage in Sublineages:
				ECDCVariants.append(VirusVariant+' ('+ECDCClassification+')')

		if query_mutations and LineageMutations:

			query_mutations = set([x.strip(' ') for x in query_mutations.split(',')])

			if not query_pango_lineage:
				query_pango_lineage = 'any_pango_lineage'

			if query_pango_lineage in LineageMutations.keys():
				required_mutations = set(LineageMutations[query_pango_lineage].split('+'))
				
				if len(required_mutations - query_mutations) == 0:
					ECDCVariants.append(VirusVariant+' ('+ECDCClassification+')')

	if VOCVOI and len(ECDCVariants) > 0:
		ECDCVariants = [ECDCVariants[-1]]

	return ','.join(ECDCVariants)


if __name__ == "__main__":

	#df=pd.read_csv(url,index_col='VirusVariant')
	df_ECDC_variants=pd.read_csv(url_static,index_col='VirusVariant')

	#DescribeVariantTable(df_ECDC_variants)

	SingleQuery(df_ECDC_variants)

	ECDCVariants = SingleQuery(df_ECDC_variants,'Q.7')
	print(ECDCVariants+'\n')

	ECDCVariants = SingleQuery(df_ECDC_variants,'HK.22')
	print(ECDCVariants+'\n')

	#pango lineage with comma separated mutations that does not meet requirement
	ECDCVariants = SingleQuery(df_ECDC_variants,'FD.1','Spike_F486P, Spike_F490S')
	print(ECDCVariants+'\n')

	#pango lineage with comma separated mutations that meets requirement
	ECDCVariants = SingleQuery(df_ECDC_variants,'FD.1','Spike_F486P, Spike_F456L, Spike_F490S')
	print(ECDCVariants+'\n')