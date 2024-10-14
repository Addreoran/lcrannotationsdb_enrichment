import requests


def read_proteins(file):
    uniprots = set()
    with open(file) as f:
        for l in f:
            if l.strip():
                uniprot = l.split()[0]
                uniprots.add(uniprot)
    return uniprots


def get_GO_quickgo(uniprots, GO_ids):
    go_quickgo = {i: set() for i in uniprots}
    for uniprot in uniprots:
        url = f"https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?includeFields=goName&includeFields=taxonName&includeFields=name&geneProductId={uniprot}"
        request = requests.get(url, headers={'Accept': 'text/tsv'})
        for line in request.text.split("\n"):
            if "GENE" not in line:
                line = line.split()
                if line:
                    if line[1] == uniprot:
                        go_quickgo[uniprot].add(line[4])
                        GO_ids.add(line[4])
    return go_quickgo


def get_annotations_GO_lcrannotdb(uniprots, go_ids):
    # info = {}
    go_lcrannotationsdb = {i: set() for i in uniprots}
    info = {i: "" for i in uniprots}
    url = f"https://lcrannotdb.lcr-lab.org/api/categories/"
    req = {
        "result_data": "annotations",
        "columns": ['UniProtACC', 'Organism', 'Annotation', 'Annotation ID', "LCR ID", "Gene ontology ID of category"],
        "search_by_column": {
            "UniProtACC": list(uniprots)
        },
    }
    request = requests.post(url, json=req)
    for l in request.text.split("\n"):
        if l.strip() and "UniProtACC" not in l:
            line = l.split(";")
            if line[0] in uniprots:
                if "GO" in line[-1]:
                    info[line[0]] = line[1]
                    go_lcrannotationsdb[line[0]].add(line[-1].strip())
                    go_ids.add(line[-1].strip())
    return go_lcrannotationsdb, info


def save_file(file, result_lcrannotdb, result_quickgo, go_family, info_organism, go_names):
    all_go_protein = set()
    all_go_categories = set()
    summary_go_proteins = {}
    summary_go_categories = {}
    summary_go_categories_no_go_protein = {}
    summary_go_protein_no_go_categories = {}
    go_family = set(go_family.split(','))
    with open(file, "w") as f:
        f.write(
            "uniprot;organism;go_protein;"
            "go_categories;"
            "go_protein-go_categories;"
            "go_categories-go_protein;"
            "common_protein_categories (union(go_protein, go_categories));"
            "go_family;go_protein-go_family;go_categories-go_family;go_categories_no_go_protein;common_protein_categories-go_family\n")
        for unirpto, goes in result_lcrannotdb.items():
            go_protein = ','.join(list(result_quickgo[unirpto]))
            go_categories = ','.join(list(goes))
            for go in goes:
                if go not in summary_go_proteins:
                    summary_go_proteins[go] = 0
                if go not in summary_go_categories:
                    summary_go_categories[go] = 0
                if go not in summary_go_categories_no_go_protein:
                    summary_go_categories_no_go_protein[go] = 0
                if go not in summary_go_protein_no_go_categories:
                    summary_go_protein_no_go_categories[go] = 0
                summary_go_categories[go] += 1
                if go not in result_quickgo[unirpto]:
                    summary_go_categories_no_go_protein[go] += 1
            for go in result_quickgo[unirpto]:
                if go not in summary_go_proteins:
                    summary_go_proteins[go] = 0
                if go not in summary_go_categories:
                    summary_go_categories[go] = 0
                if go not in summary_go_categories_no_go_protein:
                    summary_go_categories_no_go_protein[go] = 0
                if go not in summary_go_protein_no_go_categories:
                    summary_go_protein_no_go_categories[go] = 0
                if go not in goes:
                    summary_go_protein_no_go_categories[go] += 1
                summary_go_proteins[go] += 1
            common_protein_categories = ','.join(list(goes.intersection(result_quickgo[unirpto])))
            go_categories_no_go_family = ','.join(list(goes - go_family))
            common_protein_categories_no_go_family = ','.join(
                list(goes.intersection(result_quickgo[unirpto]) - go_family))
            go_protein_no_go_family = ','.join(list(result_quickgo[unirpto] - go_family))
            go_categories_no_go_protein = ','.join(list(goes - result_quickgo[unirpto]))
            go_protein_no_categories = ",".join(list(result_quickgo[unirpto] - goes))
            go_categories_no_protein = ",".join(list(goes - result_quickgo[unirpto]))

            go_family_text = ','.join(list(go_family))
            if "Accession" not in unirpto:
                f.write(f"{unirpto};"
                        f"{info_organism[unirpto]};"
                        f"{go_protein};"
                        f"{go_categories};"
                        f"{go_protein_no_categories};"
                        f"{go_categories_no_protein};"

                        f"{common_protein_categories};"
                        f"{go_family_text};"
                        f"{go_protein_no_go_family};"
                        f"{go_categories_no_go_family};"
                        f"{go_categories_no_go_protein};"

                        f"{common_protein_categories_no_go_family};"
                        f"\n")
            all_go_protein = all_go_protein.union(result_quickgo[unirpto])
            all_go_categories = all_go_categories.union(goes)
        all_common_protein_categories = ','.join(list(all_go_protein.intersection(all_go_categories)))
        go_categories_no_go_family = ','.join(list(all_go_categories - go_family))
        all_common_protein_categories_no_go_family = ','.join(
            list(all_go_protein.intersection(all_go_categories) - go_family))
        all_go_protein_no_go_family = ','.join(list(all_go_protein - go_family))
        go_categories_no_go_protein = ','.join(list(all_go_categories - all_go_protein))
        all_go_protein_no_categories = ",".join(list(all_go_protein - all_go_categories))
        all_go_categories_no_protein = ",".join(list(all_go_categories - all_go_protein))
        f.write(f"summary;;"
                f"{','.join(list(all_go_categories))};"
                f"{','.join(list(all_go_protein))};"
                f"{all_common_protein_categories};"
                f"{go_family_text};"
                f"{all_go_protein_no_categories};"
                f"{all_go_categories_no_protein};"
                f"{all_go_protein_no_go_family};"
                f"{go_categories_no_go_family};"
                f"{go_categories_no_go_protein};"
                f"{all_common_protein_categories_no_go_family};\n"
                )
        #     summary_go_proteins = {}
        #     summary_go_categories = {}
        #     summary_go_categories_no_go_protein = {}
        #     summary_go_protein_no_go_categories = {}
        f.write(
            "go_name;go_names;proteins_with_go;categories_with_go;go_categories_no_go_protein;go_protein_no_categories;\n")
        for go, no_proteins_go in summary_go_proteins.items():
            f.write(
                f"{go};{go_names[go]};{no_proteins_go};{summary_go_categories[go]};{summary_go_categories_no_go_protein[go]};{summary_go_protein_no_go_categories[go]}\n")


def get_GO_names(GO_ids):
    go_names = {i: set() for i in GO_ids}
    for go_id in GO_ids:
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"
        request = requests.get(url, headers={'Accept': 'application/json'})
        print(request.json())
        for res in request.json()["results"]:
            if res["id"] == go_id:
                go_names[go_id] = res["name"]
    return go_names
    # return go_quickgo


if __name__ == "__main__":
    rev = read_proteins("./data/rev")
    rev_family = "GO:0006355,GO:0003700,GO:0042025"
    rev_GO_ids = set(rev_family.split(","))
    rev_result_quickgo = get_GO_quickgo(rev, rev_GO_ids)
    rev_result_lcrannotationsdb, info = get_annotations_GO_lcrannotdb(rev, rev_GO_ids)
    go_names = get_GO_names(rev_GO_ids)
    save_file("./data/rev_go.csv", rev_result_quickgo, rev_result_lcrannotationsdb, rev_family, info, go_names)

    rdrp = read_proteins("./data/rdrp")
    rdrp_family = "GO:0003968,GO:0003723,GO:0039694"
    rdrp_GO_ids = set(rdrp_family.split(","))
    rdrp_result_quickgo = get_GO_quickgo(rdrp, rdrp_GO_ids)
    rdrp_result_lcrannotationsdb, info = get_annotations_GO_lcrannotdb(rdrp, rdrp_GO_ids)
    go_names = get_GO_names(rdrp_GO_ids)
    save_file("./data/rdrp_go.csv", rdrp_result_quickgo, rdrp_result_lcrannotationsdb, rdrp_family, info, go_names)

    tat = read_proteins("./data/tat")
    tat_family = "GO:0050434,GO:0001070,GO:0042025"
    tat_GO_ids = set(tat_family.split(","))
    tat_result_quickgo = get_GO_quickgo(tat, tat_GO_ids)
    tat_result_lcrannotationsdb, info = get_annotations_GO_lcrannotdb(tat, tat_GO_ids)
    go_names = get_GO_names(tat_GO_ids)
    save_file("./data/tat_go.csv", tat_result_quickgo, tat_result_lcrannotationsdb, tat_family, info, go_names)
