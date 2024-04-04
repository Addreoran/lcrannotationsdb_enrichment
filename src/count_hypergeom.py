from scipy.stats import hypergeom

import os

import requests
from django.core.management.base import BaseCommand

from lcrannotationdb.models import Protein
from statsmodels.stats.multitest import multipletests


# from scipy.stats import hypergeom


class Command(BaseCommand):
    # python manage.py initdb <seg_parameters>.
    help = "Wywołaj to do obliczeni astatystyk."

    def handle(self, *args, **options):
        folder = "./data_res/"
        new_folder = "./data_res_new/"
        files = os.listdir(folder)
        all_proteins = Protein.objects.all().count()
        new_file_data = []
        for file in files:
            file_data = read_files(f"{folder}{file}")
            new_file_data += count_hypergeom_all(file_data, all_proteins)
        new_file_data = Benjamini_Hochberg(new_file_data)
        save_to_files(new_file_data, f"./hypergeom_test_total.csv")


def Benjamini_Hochberg(data):
    significance = 0.05
    p_values = {i[12] for i in data}
    result_total = {}
    pvals_test = list(p_values)
    # todo sprawdzić czy pvals ma różne wartości
    if pvals_test:
        rest = multipletests(pvals=pvals_test, alpha=significance, method="fdr_bh")
        print(rest)

    for e, i in enumerate(pvals_test):
        result_total[i] = [rest[0][e], rest[1][e], rest[1][e]<0.05]
    result = []
    for i in data:
        result.append(i + result_total[i[12]])
    return result

def save_to_files(new_file_data, path_save):
    with open(path_save, "w") as f:
        for new_line in new_file_data:
            f.write("|".join([str(i) for i in new_line]))
            f.write("\n")

def read_files(file):
    with open(file) as f:
        for line in f:
            if line.strip():
                data = line.strip().split("|")
                if "len(" not in line:
                    yield data

def count_hypergeom_all(file, all_proteins):
    total_file = []
    for line in file:
        family_annotation_proteins = line[7]
        annotation_proteins = line[4]
        family_proteins = line[5]
        total_file.append(
            line + count_hypergeom(int(all_proteins), int(family_annotation_proteins), int(annotation_proteins),
                                   int(family_proteins)))
    return total_file


def count_hypergeom(all_proteins, family_annotation_proteins, annotation_proteins, family_proteins):
    M = all_proteins
    N = annotation_proteins
    x = family_annotation_proteins
    k = family_proteins
    stat = hypergeom.sf(x - 1, M, k, N)  # (k, M, n, N)

    return [stat, stat<0.05, M, k, N, x]