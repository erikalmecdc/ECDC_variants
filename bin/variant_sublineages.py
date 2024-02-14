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
	for VirusVariant in df.index:

		Sublineages = df.loc[VirusVariant,'included sub-lineages'].split('|')
		ECDCClassification = df.loc[VirusVariant,'ECDCClassification']
		
		LineageMutations = {}

		if not pd.isnull(df.loc[VirusVariant,'LineageMutations']):
			LineageMutations =  json.loads(df.loc[VirusVariant,'LineageMutations'])

		LineageMutationsExamples = ''
		if LineageMutations:
			i = 0

			LineageMutationsExamples += ' and additional '+str(len(LineageMutations.keys()))+' sub-lineages that fulfils the mutation criterium, e.g. '
			for lineage, mutations in LineageMutations.items():
				mutations = mutations.split('+')

				lineage.replace('any_pango_lineage','any pango lineage')

				LineageMutationsExamples+= ' '+lineage+' with '+'+'.join(mutations)+','
				#print(lineage)
				#print(mutations)
				i += 1

				if i >= limit:
					LineageMutationsExamples = LineageMutationsExamples.rstrip(',')+'.'
					break 

		print(LineageMutationsExamples)






if __name__ == "__main__":

	#df=pd.read_csv(url,index_col='VirusVariant')
	df=pd.read_csv(url_static,index_col='VirusVariant')

	DescribeVariantTable(df)




