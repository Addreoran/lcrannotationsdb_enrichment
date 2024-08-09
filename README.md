# Enrichment of LCR annotations in protein families
LCRs are usually observed as short fragments of the whole protein sequence. To analyze how functional annotations of LCRs are related to a general protein function, we analyzed enrichment of LCR annotations in InterPro families.
This method was implemented to asses usability of LCRAnnotationsDB. 

# Technologies
Project is created with Python3.6 and several packages to download data to analysis and run statistic tests (see file requirements.txt).

## Setup and run
This project was implemented in linux dystrybution and we present only setup with pip package. To run this analysis, run following steps:

1. Install Python
```
$sudo apt-get update
$sudo apt-get install python3.6
```
2. Install requirements
```
$pip install -r requirements.txt
```
3. Run analysis
```
 $python3 ./src/download_interpro_families.py --gt 70 --path "<path_to_save_tmp_files>"
 $python3 ./src/count_hypergeom.py --gt 70 --path "<path_to_save_tmp_files>"
 $python ./src/famiLy_GO_analyse.py
```
In case of our analysis, we used parameter --gt greater then 70%. It means that in this analysis, there are used only annotations that are covered with and covers LCRs in more then 70%. Parameter --path "path_to_save_tmp_files" requires path to empty folder.

File "family_GO_analyse.py" needs data folder that contains rows downloaded from InterPro database to families RdRp, REV and Tat. Input files used in our analyses are attached in data folder.
