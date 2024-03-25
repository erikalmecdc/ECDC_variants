# ECDC_variants

This repository is aimed to fascilitate ECDC SARS-CoV-2 variant classification by Pangolin lineage and/or amino-acid substitutions for internal purposes and EpiPulse/TESSy reporting.

There is one main script in Python provided with a few examples that can easily become modified as appropriate. The script either processes a single query or arbitrary GISAID EpiCoV "Patient status metadata" if placed in the "input" folder.
A resulting excel sheet will be available in the "output" folder in addition to an EpiPulse/TESSy reporting sheet.

Please note that the script reads data from a static link provided by ECDC, however this is currently (2024-03-25) only available on test data.

You may run mkdirs.sh to make the folders necessary.

Please contact Olov.Svartstrom@ecdc.europa.eu in case there are questions/issues.
