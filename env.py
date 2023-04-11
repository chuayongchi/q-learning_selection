import numpy as np

class GeneEnv():

  def __init__(self, ggc_df, top_gene_id_list, count_row, gene_id, sample):
    self.gene_num = count_row
    self.gene_list = gene_id

    self.terminated_factor = False
    self.selection = np.zeros(shape=(self.gene_num, ), dtype=int)

    self.top_gene_list = top_gene_id_list
    self.ggc = ggc_df.values.tolist()
    self.index_top_genes = np.zeros(shape=(10, ), dtype=int)
    self.expression = sample

  def is_connective(self):
    is_connective = True
    n = len(self.selection)
    
    for i in range (n):
      for j in range (n):
        if ((self.selection[i] == 1) and (self.selection[j] == 1) and (i != j)):
          if (self.ggc[i][j] < 0.1):
            is_connective = False

    print("is_connective: ", is_connective)
    return is_connective

  def play_step(self, action, terminated_index):
    expressed = 0

    if (action == terminated_index):
      self.terminated_factor = True
      expressed = 0

    else:
      self.selection[action] = 1
      expressed = self.expression[action]

    return self.terminated_factor, self.selection, expressed

  def reset(self):
    self.terminated_factor = False
    self.selection = np.zeros(shape=(self.gene_num, ), dtype=int)
  
  def get_index_top_genes(self):
    num_top_gene_list = len(self.top_gene_list)
    print("num_top_gene_list: ", num_top_gene_list)
    for i in range (num_top_gene_list):
      for j in range (self.gene_num):
        if(self.top_gene_list[i] == self.gene_list[j]):
          self.index_top_genes[i] = j
    return self.index_top_genes

  def cal_max_connectivity(self, selection):
    connectivity = 0
    n = len(selection)
    hub_gene = []
    temp = []
    for i in range (n):
      for j in range (n):
        if ((selection[i] == 1) and (selection[j] == 1) and (i != j)):
          connectivity = connectivity + self.ggc[i][j]
        elif ((selection[i] == 1) and (selection[j] == 1) and (i == j)):
          connectivity = connectivity + 0
      hub_gene.append(selection[i])
      temp.append(connectivity)
      connectivity = 0
    
    connectivity = temp[np.argmax(temp)]
    return connectivity