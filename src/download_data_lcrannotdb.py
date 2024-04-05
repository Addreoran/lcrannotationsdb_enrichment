class LCRAnnotDBData:
    def __init__(self, path, gt):
        self.path = path
        self.protein_no = None
        self.gt=gt

    def get_categories_pk_list(self):
        pass

    def get_annotations_category(self, category_pk):
        # {annotation:{"uniprotID":"", SourceID:"", SourceName:"", Name:"", pk:""}}
        pass

    def get_category_name_by_pk(self, pk):
        pass

    def get_proteins_by_uniprotAccs(self, Uniprots):
        # {"Uniprot":"", "Header":"", "Sequence":""}
        pass

