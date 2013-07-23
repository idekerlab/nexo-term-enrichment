# Term Enrichement
import numpy as np
import scipy as sp
from scipy import stats

from flask import Flask
from flask.ext import restful

app = Flask(__name__)
api = restful.Api(app)

class TermMapping:

  __gene2nexo_terms = {}
  __nexo_term2genes = {}

  def __init__(self):
    print "Loading resource"
    self.__gene2nexo_terms = self.createMap("data/gene2terms.txt")
    self.__nexo_term2genes = self.createMap("data/term2genes.txt")

  def createMap(self, file_name):
    mapping = {}
    with open(file_name, "r") as f:
      for line in f:
        line = line.rstrip()
        parts = line.split("\t")
        geneID = parts[0]
        terms = parts[1].split(",")
        mapping[geneID] = terms
    f.close()
    return mapping


  def getMapping(self):
    return self.__gene2nexo_terms

  def getTermMapping(self):
    return self.__nexo_term2genes


class HypergeometricTester:
  def __init__(self):
    print("Tester Initialized")

  def performTest(self,gene2terms,term2genes,genes_of_interest,total_num_genes):
    sample_map = self.__calculateTermFrequency(genes_of_interest, gene2terms)

    p_values = {}

    n = len(genes_of_interest)

    for term in sample_map:
      # Calculate p-value
      sampled = sample_map[term]
      assigned_genes = term2genes[term]
      k = len(sampled)
      m = len(assigned_genes)
      p = stats.hypergeom.pmf(k,total_num_genes,m,n)
  

  def __calculateTermFrequency(self, genes_of_interest, gene2terms):
    sample_term_map = {}

    for gene in genes_of_interest:
      terms = gene2terms[gene]
      for term in terms:
        if term in sample_term_map:
          assigned_genes = sample_term_map[term]
          assigned_genes.add(gene)
          sample_term_map[term] = assigned_genes
        else:
          assigned_genes = set([])
          assigned_genes.add(gene)
          sample_term_map[term] = assigned_genes
    
    return sample_term_map




class TermEnrichement(restful.Resource):
  def get(self):
    p = stats.hypergeom.pmf(3,6381,27,5)
    return {'nexo': 'v1', 'P-val': p}

api.add_resource(TermEnrichement, '/')

# Main 
if __name__ == '__main__':
  termMapper = TermMapping()
  print(termMapper.getTermMapping())
  app.run(debug=True)

