# Table of Contents

- [Description](#description)

    - [Enrichment of LCR annotations in protein families](#enrichment-of-lcr-annotations-in-protein-families)
    - [Gene Ontology terms of LCR annotation categories in selected protein families](#gene-ontology-terms-of-lcr-annotation-categories-in-selected-protein-families)

- [Requirements](#requirements-)

- [Setup and run](#requirements-)


# Description

This repository contains code to generate data for publication "LCRAnnotationsDB: a database of low complexity regions
functional and structural annotations":

1. Appendix A - it contains data about enrichment of LCR annotations in protein families (see section 2.3 of manuscript),
2. Appendix D - it contains data with Gene Ontologies of selected families with proteins enriched in low complexity
   regions, these data are detailed described in section 3.1.1 of manuscript.

## Enrichment of LCR annotations in protein families

First, the user needs to download information about InterPro Families from the InterPro database. 
To get it, it is necessary to run the script "./src/download_interpro_families.py". It needs parameters "--gt" to 
to obtain families of only proteins with categories with annotations covering all LCRs with value greater 
than parameter "--gt" and vice versa. In the case of the analyses in the manuscript, we used the value 70 of the parameter "--gt".
The output folder "./data_res/" contains detailed information about overlap between proteins with annotations 
and protein families. 

 In the next step, the script "./src/count_hypergeom.py" calculates the hypergeometric test with Benjamini-Hochberg correction. 
 This script requires files from the "./data_res/" folder from the previous script. In this analysis we have excluded annotations with 
 with InterPro Family IDs from the LCR dataset. The output from this file is saved in the file "./hypergeom_test_total.csv".


```
 $python3 ./src/download_interpro_families.py --gt 70
 $python3 ./src/count_hypergeom.py --gt 70
```


The scripts included in this repository contains workflow of the analysis presented here due to the considerable volume 
of data required. To obtain the results of this analysis, we utilised other scripts with a direct connection to the 
database on the server. Due to the lack of access to this server for users, we rewrite it in order to 
visualise the execution of this analysis. With very long time of code execution, users can get correct result, 
with the new version of InterPro Database.

## Gene Ontology terms of LCR annotation categories in selected protein families

As an output, there is a summary of the Gene Ontology (GO) terms associated with the families, proteins, and categories 
analysed in the use case of manuscript about LCRAnnotationsDB (section 3.1.1). The GOs originate from the QuickGO database, while those pertaining to the categories have
been sourced from the LCRAnnotationsDB database. The information concerning the families has been obtained from the 
InterPro database.

The data were collated and formatted into Excel sheets, titled "REV - GO of families, categories and proteins," 
"RdRp - GO of families, categories and proteins," and "Tat - GO of families, categories and proteins." 
These sheets can be found in Appendix D of the publication 
"LCRAnnotationsDB: a database of low complexity regions functional and structural annotations" 
from which they were also used to create Appendix D.

After the setup commands, this part can be run with command:

```
 $python3 ./src/famiLy_GO_analyse.py
```

This script downloads Gene Ontology (GO) terms for proteins from the QuickGO database and GO terms categories of LCR 
annotations from LCRAnnotationsDB. The proteins used in this analysis are sourced from the folder "./Data/". Finally,
the coverage of families, GO terms and GO terms of proteins and GO terms of categories of LCR annotations for each 
protein are compared. 

# Requirements

Project is created with Python3.6 and several packages to download data to analysis and run statistic tests (see file
requirements.txt).

# Setup and run

This project was implemented in the Linux distribution and we present only the setup with the pip package. To perform this analysis
the following steps:

1. Install Python

```
$sudo apt-get update
$sudo apt-get install python3.6
```

2. Install requirements

```
$pip install -r requirements.txt
```

3. Run analysis - shortcut

Here you can see general commands to execute scripts from this project. To see a more detailed description of the execution analysis, see the [Description](#description) chapter.

```
 $python3 ./src/download_interpro_families.py --gt 70"
 $python3 ./src/count_hypergeom.py --gt 70"
 $python3 ./src/famiLy_GO_analyse.py
```

In the case of our analysis, we used the parameter --gt greater than 70%. This means that in this analysis only
annotations that are more than 70% covered by LCRs. 

The file './src/family_GO_analyse.py' requires a data folder containing rows from the InterPro database for the
RdRp, REV, and Tat families. The input files used in our analyses are included in the data folder.



