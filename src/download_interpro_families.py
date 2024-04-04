import os

import requests
from django.core.management.base import BaseCommand

from lcrannotationdb.models import DescriptionLCRRelation, CategoriesDescriptionRelation, Categories, \
    Protein

class Command(BaseCommand):
    help = "Wywołaj to do obliczeni astatystyk."

    def add_arguments(self, parser):
        parser.add_argument(
            '--gt',
            action='store',
            default=None,
            help='Delete poll instead of closing it',
        )

        parser.add_argument(
            '--file',
            action='store',
            default=None,
            help='Delete poll instead of closing it',
        )


    def handle(self, *args, **options):
        self.create_table(options)

    def create_table(self, options):
        proteins_no = Protein.objects.all().count()
        test_multi = {}

        annotations = DescriptionLCRRelation.objects.filter(
            LCRCoverage__gte=float(options["gt"]),
            DescriptionCoverage__gte=float(options["gt"])
        ).values_list("DescriptionID", flat=True)
        categories = CategoriesDescriptionRelation.objects.filter(DescriptionID__in=annotations)  # .filter(
        # CategoriesID__Name__icontains="rna binding")

        category_once = Categories.objects.filter(
            CategoriesID__in=categories.values_list("CategoriesID", flat=True)).distinct()
        biggest_cat = {}
        nr_cat = category_once.count()
        for e, category_pk in enumerate(category_once):
            print(e, nr_cat)
            annotations = DescriptionLCRRelation.objects.filter(
                LCRCoverage__gte=float(options["gt"]),
                DescriptionCoverage__gte=float(options["gt"])
            ).filter(DescriptionID__categoriesdescriptionrelation__CategoriesID__pk=category_pk.pk)
            biggest_cat[category_pk] = annotations.count()
        # with open(options["file"], "w") as f:
        tmp_families_per_protein = self.read_tmp_data("tmp_families_per_protein.csv")
        tmp_families_proteins = self.read_tmp_data("tmp_families_proteins.csv")
        meaning_result = []
        for category_pk in sorted(biggest_cat.items(), key=lambda x: x[1], reverse=True):
            all_data = {}
            category_pk = category_pk[0]
            category = CategoriesDescriptionRelation.objects.filter(CategoriesID__pk=category_pk.pk)
            annotations = DescriptionLCRRelation.objects.filter(
                LCRCoverage__gte=float(options["gt"]),
                DescriptionCoverage__gte=float(options["gt"])
            ).filter(DescriptionID__categoriesdescriptionrelation__CategoriesID__pk=category_pk.pk)
            print("annotations")
            ann_nr = annotations.count()
            for e, annotation in enumerate(annotations):
                print(e, ann_nr)
                if annotation.LcrID.ProteinID.UniprotID:
                    if annotation.LcrID.ProteinID.UniprotID not in all_data:
                        all_data[annotation.LcrID.ProteinID.UniprotID] = {"annotations": set(), "family": set()}
                    if annotation.DescriptionID.SourceID is None:
                        all_data[annotation.LcrID.ProteinID.UniprotID]["annotations"].add(
                            f"{annotation.DescriptionID.SourceName.SourceName}={annotation.DescriptionID.Description}")
                    else:
                        if annotation.DescriptionID.SourceID != annotation.LcrID.ProteinID.UniprotID:
                            all_data[annotation.LcrID.ProteinID.UniprotID]["annotations"].add(
                                annotation.DescriptionID.SourceID)
                        else:
                            all_data[annotation.LcrID.ProteinID.UniprotID]["annotations"].add(
                                f"{annotation.DescriptionID.SourceName.SourceName}={annotation.DescriptionID.Description}")
                else:
                    print("brak Uniprotid", annotation.LcrID.ProteinID.pk)
            nr_all_data = len(all_data)
            counter_all_data = 0
            tmp_proteins = []
            for protein, data_protein in all_data.items():
                counter_all_data += 1
                if protein not in tmp_families_per_protein:
                    tmp_families_per_protein[protein] = self.get_interpro_family(protein)
                    all_data[protein]["family"] = tmp_families_per_protein[protein]
                else:
                    all_data[protein]["family"] = tmp_families_per_protein[protein]
            families_proteins = {}
            counter_all_data = 0
            for prot, prot_data in all_data.items():
                counter_all_data += 1
                for family in all_data[prot]["family"]:
                    if family not in tmp_families_proteins:
                        families_proteins[family] = self.analyse_interpro_family(family)
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
                        tmp_families_proteins[annotation] = self.analyse_interpro_family(annotation)
                        annotations_proteins_interpro[annotation] = tmp_families_proteins[annotation]
                    else:
                        annotations_proteins_interpro[annotation] = tmp_families_proteins[annotation]
            res = []
            for annotation, annot_proteins in annotations_proteins.items():
                for family, family_proteins in families_proteins.items():
                    if len(annot_proteins & family_proteins) > 0:
                        res.append([str(category_pk.pk),
                                    str(Categories.objects.get(pk=category_pk.pk).Name),
                                    str(annotation),
                                    str(family),
                                    str(len(annot_proteins)),
                                    str(len(family_proteins)),
                                    str(len(annotations_proteins_interpro.get(annotation, []))),
                                    str(len(annot_proteins & family_proteins)),
                                    str(len(annot_proteins & family_proteins) / len(annot_proteins)),
                                    str(len(annot_proteins & family_proteins) / len(family_proteins)),
                                    str(len(annotations_proteins_interpro.get(annotation, [])) / len(annot_proteins)),
                                    str(len(annotations_proteins_interpro.get(annotation, [])) / len(family_proteins)),
                                    ])
                        if (len(annot_proteins & family_proteins) / len(annot_proteins)) > 0.20 and (len(
                                annot_proteins & family_proteins) / len(
                            family_proteins)) > 0.20 and len(
                            annot_proteins) > 10 and len(family_proteins) > 10:
                            meaning_result.append([str(category_pk.pk),
                                                   str(Categories.objects.get(pk=category_pk.pk).Name),
                                                   str(annotation),
                                                   str(family),
                                                   str(len(annot_proteins)),
                                                   str(len(family_proteins)),
                                                   str(len(annotations_proteins_interpro.get(annotation, []))),
                                                   str(len(annot_proteins & family_proteins)),
                                                   str(len(annot_proteins & family_proteins) / len(
                                                       annot_proteins)).replace(".", ","),
                                                   str(len(annot_proteins & family_proteins) / len(
                                                       family_proteins)).replace(".", ","),
                                                   str(len(annotations_proteins_interpro.get(annotation, [])) / len(
                                                       annot_proteins)).replace(".", ","),
                                                   str(len(annotations_proteins_interpro.get(annotation, [])) / len(
                                                       family_proteins)).replace(".", ","),
                                                   ])
                            self.save_as_fasta(path=f"./data_res_proteins/an_{annotation.replace('/', '_')}",
                                               dataset=self.select_annotations_sequences(annot_proteins))
                            self.save_as_fasta(path=f"./data_res_proteins/fam_{annotation.replace('/', '_')}",
                                               dataset=self.select_family_sequences(family_proteins))
                            self.save_as_fasta(path=f"./data_res_proteins/no_an_fam_{annotation.replace('/', '_')}",
                                               dataset=self.select_family_sequences(
                                                   family_proteins - annot_proteins))
                            self.save_as_fasta(path=f"./data_res_proteins/an_no_fam_{annotation.replace('/', '_')}",
                                               dataset=self.select_family_sequences(
                                                   annot_proteins - family_proteins))
            with open(f"./data_res/{str(Categories.objects.get(pk=category_pk.pk).Name).replace('/', '_')}", "w") as f:
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
        with open(f"./result_general.csv", "w") as f:
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
            if meaning_result:
                for i in meaning_result:
                    if i:
                        f.write("|".join(i))
                        f.write("\n")
                        print("geeral result", "|".join(i))


    def save_as_fasta(self, path, dataset):
        with open(path, "w") as f:
            for header, header_sequence in dataset.items():
                header = header_sequence[0]
                sequence = header_sequence[1]
                f.write(f"{header}\n")
                f.write(f"{sequence}\n")

    def select_family_sequences(self, UniprotAccs):
        fasta_database = {}
        proteins_db = Protein.objects.filter(UniprotID__in=UniprotAccs)
        no_in_db = UniprotAccs - set(proteins_db.values_list("UniprotID", flat=True))
        for protein in proteins_db:
            fasta_database[protein.UniprotID] = (protein.Header, protein.Sequence)
        for protein_noid in no_in_db:
            fasta_database[protein_noid] = self.get_sequence(protein_noid)
        return fasta_database

    def select_annotations_sequences(self, UniprotAccs):
        fasta_database = {}
        proteins_db = Protein.objects.filter(UniprotID__in=UniprotAccs)
        for protein in proteins_db:
            fasta_database[protein.UniprotID] = (protein.Header, protein.Sequence)
        return fasta_database

    def get_sequence(self, uniprot):
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta"
        req = requests.get(url)
        if req.text:
            text_req = req.text.split("\n")
            header = text_req[0]
            res = text_req[1:]
        return header, "/n".join(res)

    def analyse_interpro_family(self, interpro_family):
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

    def read_tmp_data(self, file):
        result = {}
        if os.path.exists(file):
            with open(file) as f:
                for i in f:
                    if i.strip():
                        key_item, value_item = i.strip().split("|")
                        result[key_item] = set([i for i in value_item.split(";") if i])
        return result

def get_interpro_family(self, protein):
        family_proteins = set()
        nr = 0
        # for family in interpro_families:
        request_url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{protein}?format=json"
        while request_url:
            try:
                print(request_url)
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

    return interpro_family, family_proteins
