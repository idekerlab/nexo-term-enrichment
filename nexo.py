# Term Enrichement
import numpy as np
import scipy as sp
from scipy import stats

from flask import Flask
from flask.ext import restful

app = Flask(__name__)
api = restful.Api(app)

class TermMapping:

  __gene2nexo_term = {}

  def __init__(self):
    print "Loading resource"
    with open("data/gene2terms.txt", "r") as f:
      for line in f:
        line = line.rstrip()
        parts = line.split("\t")
        geneID = parts[0]
        terms = parts[1].split(",")
        self.__gene2nexo_term[geneID] = terms

  def getMapping(self):
    return self.__gene2nexo_term

class HypergeometricTester:
  def __init__(self):
    print("Tester Initialized")

  def performTest(self, genes_of_interest):
    num_genes = len(genes_of_interest)


class TermEnrichement(restful.Resource):
  def get(self):
    return {'nexo': 'v1'}

api.add_resource(TermEnrichement, '/')

# Main 
if __name__ == '__main__':
  termMapper = TermMapping()
  app.run(debug=True)

