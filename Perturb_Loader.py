
import anndata
import gdown


def download_perturb_data():
  files = '''ann_guide<>https://drive.google.com/file/d/1f8IubR_kRwhx3js2iib2BEQxkYDDZa7t/view?usp=sharing
  ann_exp<>https://drive.google.com/file/d/1FzYOdI1iEs4z7wYe8PYOxliHeJt90AJ_/view?usp=sharing
  ann_Element_x_tested_genes<>https://drive.google.com/file/d/1--2k8r9XQFi65qzcP4fZxFX3NBNIV1tw/view?usp=sharing
  ann_Element_guide<>https://drive.google.com/file/d/1xQ9InZE1yDxjrR5HwHJVvKONel55Zaxs/view?usp=sharing'''.replace(' ', '').split('\n')
  
  for f in files:
    ann_file_name, ann_path = f .split('<>')
    ann_path = ann_path.replace('/view?usp=sharing', '').split('/')[-1]
    ann_file_name = ann_file_name + '.h5ad'
    print (ann_file_name, ann_path , 'downloading')
    gdown.download(id =ann_path, output=ann_file_name, quiet=False)




def load_perturb():

  ann_exp = anndata.read_h5ad('ann_exp.h5ad')
  ann_guide = anndata.read_h5ad('ann_guide.h5ad')
  ann_Element_guide= anndata.read_h5ad('ann_Element_guide.h5ad')
  ann_Element_x_tested_genes= anndata.read_h5ad('ann_Element_x_tested_genes.h5ad')


  p = PERTURB_MANIPULATE(ann_exp, ann_guide, ann_Element_guide, ann_Element_x_tested_genes)


  return p




class PERTURB_MANIPULATE():
  def __init__(self, cell_x_gene_ann, cell_x_guide, element_x_guide, element_x_tested_gene):
    self.cell_x_gene_ann = cell_x_gene_ann
    self.cell_x_guide = cell_x_guide
    self.element_x_guide = element_x_guide
    self.element_x_tested_gene = element_x_tested_gene
    self.check_integrity()
  
  

  def element_description(self):
    return set(self.element_x_tested_gene.obs['GUIDE_TYPE'].values)


  def extract_element_per_type(self, guide_type):
    bool_table_ele =  (self.element_x_tested_gene.obs['GUIDE_TYPE'].values == guide_type)
    return self.element_x_tested_gene.obs.index.values [bool_table_ele]


  def show_barcodes(self):
    print (f'Number of barcodes :{self.cell_x_guide.obs.index.values.shape[0]}')
    return  self.cell_x_guide.obs.index.values.tolist()
  def show_guides(self):
    print (f'Number of Guides :{self.cell_x_guide.var.index.values.shape[0]}')
    return  self.cell_x_guide.var.index.values.tolist()

  def show_elements(self):
    print (f'Number of Elements :{self.element_x_guide.obs.index.values.shape[0]}')

    return  self.element_x_guide.obs.index.values.tolist()

  def check_integrity(self):
    assert set(self.cell_x_gene_ann.obs.index.values ) == set(self.cell_x_guide.obs.index.values), 'cell_x_gene_ann and  cell_x_guide has different obs len'
 
  def capture_gene_count(self, gene):
    return self.cell_x_gene_ann[:, gene].to_df().values

  def capture_guide(self, guide_name):   
    return self.cell_x_guide[:, guide_name].to_df().values.reshape(-1)

  def capture_element_guides(self, element_name):  
    return  self.element_x_guide.var.index.values[(self.element_x_guide[element_name:, ].to_df().values == 1)[0]].tolist()

  def capture_covariates(self):
    return self.cell_x_guide.obs

  def capture_tested_genes_vector(self, guide_name):
    return self.element_x_tested_gene[guide_name,:].to_df().values[0]

  def capture_tested_genes_names(self, guide_name):
    genes_names =  self.element_x_tested_gene.var.index.values
    gene_presence =  self.capture_tested_genes_vector(guide_name)
    return   genes_names[gene_presence == 1]

  def genes_to_test_count(self, guide_name):
    return self.cell_x_gene_ann[:, self.capture_tested_genes_names(guide_name)].to_df()

  #p = PERTURB_MANIPULATE(ann_exp, ann_guide, ann_Element_guide, ann_Element_x_tested_genes)
