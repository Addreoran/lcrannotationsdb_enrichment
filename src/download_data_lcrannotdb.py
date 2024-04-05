import requests


class LCRAnnotDBData:
    def __init__(self, path, gt):
        self.path = path
        self.protein_no = None
        if self.protein_no is None:
            self.protein_no = self.get_protein_no()
        self.gt = gt
        self.categories_data = None

    def get_protein_no(self):
        url = f"https://lcrannotdb.lcr-lab.org/api/proteins/"
        req = {
            "result_data": "proteins",
            "columns": ['UniProtACC'],
        }
        request = requests.post(url, json=req)
        result = set()
        for line in request.text.split("\n"):
            if "UniProtACC" not in line:
                line = line.split(";")
                result.add(line[0])
        return len(result)

    def get_categories(self):
        url = f"https://lcrannotdb.lcr-lab.org/api/categories/"
        req = {
            "result_data": "categories",
            "columns": ['Category ID', 'Category'],
        }
        request = requests.post(url, json=req)
        for line in request.text.split("\n"):
            if "Category ID;Category" not in line:
                line = line.split(";")
                self.categories_data[int(line[0])] = line[1]

    def get_categories_pk_list(self):
        if self.categories_data is None:
            self.get_categories()
        return list(self.categories_data.keys())

    def get_category_name_by_pk(self, pk):
        return self.categories_data[int(pk)]

    def get_proteins_by_uniprotAccs(self, Uniprots):
        # [{"Uniprot":"", "Header":"", "Sequence":""}]
        url = f"https://lcrannotdb.lcr-lab.org/api/proteins/"
        req = {
            "result_data": "annotations",
            "columns": ['UniProtACC', 'Protein header', 'Protein Sequence'],
            "search_by_column": {
                "UniProtACC": Uniprots
            },
        }
        request = requests.post(url, json=req)
        result = []
        for line in request.text.split("\n"):
            if "UniProtACC;Protein header;Protein Sequence" not in line:
                line = line.split(";")
                result.append(
                    {"Uniprot": line[0], "Header": line[1], "Sequence": line[2]}
                )
        return result

    def get_annotations_category(self, category_pk):
        # [{"uniprotID":"", SourceID:"", SourceName:"", Name:"", pk:""}]
        url = f"https://lcrannotdb.lcr-lab.org/api/categories/"
        req = {
            "result_data": "annotations",
            # "LCRCoverage":self.gt,
            # "DescriptionCoverage":self.gt,
            "columns": [
                # 'UniProtACC', 'Source ID', 'Source database', 'Annotation',
                "Annotation ID"
            ],
            "search_by_column": {
                "Category ID": [category_pk]
            },
        }

        request = requests.post(url, json=req)
        annotations_ids = set()
        for line in request.text.split("\n"):
            if "Annotation ID" not in line:
                line = line.split(";")
                annotations_ids.add(
                    line[0]
                )
        # return result
        url = f"https://lcrannotdb.lcr-lab.org/api/annotations/"
        req = {
            "result_data": "annotations",

            "columns": [
                'UniProtACC', 'Source ID', 'Source database', 'Annotation',
                "Annotation ID"
            ],
            "search_by_column": {
                "Annotation ID": list(annotations_ids),
                "LCRCoverage": self.gt,
                "DescriptionCoverage": self.gt,
            },
        }

        request = requests.post(url, json=req)
        result = []
        for line in request.text.split("\n"):
            if "UniProtACC;Source ID;Source database;Annotation;Annotation ID" not in line:
                line = line.split(";")
                result.append(
                    {"uniprotID": line[0], "SourceID": line[1], "SourceName": line[2], "Name": line[3], "pk": line[4]}
                )
        return result
