import os

import requests

import click

from download_data_lcrannotdb import LCRAnnotDBData


@click.command()
@click.option('--gt', default=70, help='The lowest value of LCR and annotation coverage.')
@click.option('--path', default="", help='Path to folder with data about annotations.')
def main(gt, path):
    lcrannotdb_data = LCRAnnotDBData(gt)
    category_once = lcrannotdb_data.get_categories_pk_list()
    biggest_cat = {}
    for e, category_pk in enumerate(category_once):
        annotations = lcrannotdb_data.get_annotations_category(category_pk)
        biggest_cat[category_pk] = len(annotations)
    tmp_families_per_protein = read_tmp_data("tmp_families_per_protein.csv")
    tmp_families_proteins = read_tmp_data("tmp_families_proteins.csv")
    meaning_result = []
    if not os.path.exists("./data_res_proteins/"):
        os.mkdir("./data_res_proteins/")
    if not os.path.exists("./data_res/"):
        os.mkdir("./data_res/")
    for category_pk in sorted(biggest_cat.items(), key=lambda x: x[1], reverse=True):
        all_data = {}
        category_pk = category_pk[0]
        annotations = lcrannotdb_data.get_annotations_category(category_pk)
        print("annotations")
        ann_nr = len(annotations)
        for e, annotation in enumerate(annotations):
            print(e, ann_nr)
            if annotation["UniprotID"]:
                if annotation["UniprotID"] not in all_data:
                    all_data[annotation["UniprotID"]] = {"annotations": set(), "family": set()}
                if annotation["SourceID"] is None:
                    all_data[annotation["UniprotID"]]["annotations"].add(
                        f"{annotation['SourceName']}={annotation['Name']}")
                else:
                    if annotation["SourceID"] != annotation["UniprotID"]:
                        all_data[annotation["UniprotID"]]["annotations"].add(
                            annotation["SourceID"])
                    else:
                        all_data[annotation["UniprotID"]]["annotations"].add(
                            f"{annotation['SourceName']}={annotation['Name']}")
            else:
                print("brak Uniprotid", annotation["pk"])
        nr_all_data = len(all_data)
        counter_all_data = 0
        for protein, data_protein in all_data.items():
            counter_all_data += 1
            if protein not in tmp_families_per_protein:
                tmp_families_per_protein[protein] = get_interpro_family(protein)
                all_data[protein]["family"] = tmp_families_per_protein[protein]
            else:
                all_data[protein]["family"] = tmp_families_per_protein[protein]
        families_proteins = {}
        counter_all_data = 0
        for prot, prot_data in all_data.items():
            counter_all_data += 1
            for family in all_data[prot]["family"]:
                if family not in tmp_families_proteins:
                    families_proteins[family] = analyse_interpro_family(family)
                    tmp_families_proteins[family] = families_proteins[family]
                else:
                    families_proteins[family] = tmp_families_proteins[family]
        annotations_proteins = {}
        nr_merge = 0
        for protein, data_protein in all_data.items():
            print(nr_merge, nr_all_data)
            nr_merge += 1
            for annotation in all_data[protein]["annotations"]:
                if annotation not in families_proteins:
                    if annotation not in annotations_proteins:
                        annotations_proteins[annotation] = set()
                    annotations_proteins[annotation].add(protein)
        annotations_proteins_interpro = {}
        nr_annot_len = len(annotations_proteins)
        for e, annotation in enumerate(annotations_proteins.keys()):
            if "IPR" in annotation:
                print(nr_annot_len, e, annotation)
                if annotation not in tmp_families_proteins:
                    tmp_families_proteins[annotation] = analyse_interpro_family(annotation)
                    annotations_proteins_interpro[annotation] = tmp_families_proteins[annotation]
                else:
                    annotations_proteins_interpro[annotation] = tmp_families_proteins[annotation]
        res = []

        with open(f"./data_res/{str(lcrannotdb_data.get_category_name_by_pk(category_pk)).replace('/', '_')}",
                  "w") as f:
            f.write(f"category_pk.pk|"
                    f"Categories.objects.get(pk=category_pk.pk).Name|"
                    f"annotation|"
                    f"family|"
                    f"len(annot_proteins)|"
                    f"len(family_proteins)|"
                    f"len(annotations_proteins_interpro.get(annotation, []))|"
                    f"len(annot_proteins & family_proteins)|"
                    f"len(annot_proteins & family_proteins) / len(annot_proteins)|"
                    f"len(annot_proteins & family_proteins) / len(family_proteins)|"
                    f"len(annotations_proteins_interpro.get(annotation, [])) / len(annot_proteins)|"
                    f"len(annotations_proteins_interpro.get(annotation, [])) / len(family_proteins))\n")
            if res:
                for i in res:
                    if i:
                        f.write("|".join(i))
                        f.write("\n")
                        print("|".join(i))



def save_as_fasta(path, dataset):
    with open(path, "w") as f:
        for header, header_sequence in dataset.items():
            header = header_sequence[0]
            sequence = header_sequence[1]
            f.write(f"{header}\n")
            f.write(f"{sequence}\n")


def select_family_sequences(UniprotAccs, lcrannotdb_data):
    fasta_database = {}
    proteins_db = lcrannotdb_data.get_proteins_by_uniprotAccs(UniprotAccs)
    no_in_db = UniprotAccs - set([protein["UniprotID"] for protein in proteins_db])
    for protein in proteins_db:
        fasta_database[protein["UniprotID"]] = (protein["Header"], protein["Sequence"])
    for protein_noid in no_in_db:
        fasta_database[protein_noid] = get_sequence(protein_noid)
    return fasta_database


def select_annotations_sequences(UniprotAccs, lcrannotdb_data):
    fasta_database = {}
    proteins_db = lcrannotdb_data.get_proteins_by_uniprotAccs(UniprotAccs)
    for protein in proteins_db:
        fasta_database[protein["UniprotID"]] = (protein["Header"], protein["Sequence"])
    return fasta_database


def get_sequence(uniprot):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta"
    req = requests.get(url)
    if req.text:
        text_req = req.text.split("\n")
        header = text_req[0]
        res = text_req[1:]
    return header, "/n".join(res)


def analyse_interpro_family(interpro_family):
    family_proteins = set()
    request_url = f"https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/interpro/{interpro_family}?format=json"
    while request_url:
        print(request_url)
        try:
            req = requests.get(
                request_url)
        except Exception as e:
            print(e)
            request_url = request_url
            continue
        request_url = None

        if req.status_code == 200:
            json_req = req.json()
            if json_req["next"] is not None:
                request_url = json_req["next"]
            for protein in json_req["results"]:
                family_proteins.add(protein["metadata"]["accession"])
    with open("tmp_families_proteins.csv", "a") as f:
        f.write(f"{interpro_family}|")
        f.write(";".join(family_proteins))
        f.write("\n")
    return family_proteins


def read_tmp_data(file):
    result = {}
    if os.path.exists(file):
        with open(file) as f:
            for i in f:
                if i.strip():
                    key_item, value_item = i.strip().split("|")
                    result[key_item] = set([i for i in value_item.split(";") if i])
    return result


def get_interpro_family(protein):
    family_proteins = set()
    nr = 0
    request_url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{protein}?format=json"
    while request_url:
        try:
            req = requests.get(
                request_url)
        except Exception as e:
            print(e)
            request_url = request_url
            continue
        request_url = None

        if req.status_code == 200:
            json_req = req.json()
            if json_req["next"] is not None:
                request_url = json_req["next"]
            for protein_data in json_req["results"]:
                if protein_data["metadata"]["type"] == "family":
                    family_proteins.add(protein_data["metadata"]["accession"])
    with open("tmp_families_per_protein.csv", "a") as f:
        f.write(f"{protein}|")
        f.write(";".join(family_proteins))
        f.write("\n")
    return family_proteins


if __name__ == '__main__':
    main()
