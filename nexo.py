# Term Enrichement
import numpy as np
import scipy as sp

import qvalue

from scipy import stats

from flask import Flask
from flask.ext import restful
from flask.ext.restful import Resource, fields, reqparse

import json

# Constants
GENES = 'genes'
ALPHA = 'alpha'

# p-value cutoff
p_val_threshold = 0.01

# At least this number of genes should be assigned
gene_threshold = 2

app = Flask(__name__)
api = restful.Api(app)

parser = reqparse.RequestParser()
parser.add_argument(GENES, type=str)
parser.add_argument(ALPHA, type=float)


class TermMapping:

    _gene2nexo_terms = {}
    _nexo_term2genes = {}
    _num_nexo_genes = 0

    def __init__(self):
        print "Loading resource"
        self._gene2nexo_terms = self.createMap("data/gene2terms.txt")
        self._nexo_term2genes = self.createMap("data/term2genes.txt")
        # Count genes.  Root term has all genes.
        all_genes = self._nexo_term2genes["joining_root"]
        self._num_nexo_genes = len(all_genes)
        print("NeXO Genes = " + str(self._num_nexo_genes))

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


    def get_gene_mapping(self):
        return self._gene2nexo_terms

    def get_term_mapping(self):
        return self._nexo_term2genes

    def get_gene_count(self):
        return self._num_nexo_genes


class HypergeometricTest:

    _mapper = None

    def __init__(self):
        self._mapper = TermMapping()

    def performTest(self, genes_of_interest, cutoff):
        if cutoff is None:
            cutoff = p_val_threshold

        total_num_genes = self._mapper.get_gene_count()
        term2genes = self._mapper.get_term_mapping()
        gene2terms = self._mapper.get_gene_mapping()

        sample_map = self.__calculateTermFrequency(genes_of_interest, gene2terms)

        n = len(genes_of_interest)

        # Number of tests performed: will be used for correction.
        num_tests = len(sample_map)

        results = [ {} ] * num_tests

        pvals = np.ones(num_tests)

        idx = 0
        for term in sample_map:
            # Calculate p-value
            sampled = sample_map[term]
            assigned_genes = term2genes[term]
            k = len(sampled)
            m = len(assigned_genes)

            p = stats.hypergeom.pmf(k,total_num_genes,m,n)
            p_corrected = p * num_tests
            pvals[idx] = p
            
            result = {}
            result['id'] = term
            result['p-value'] = p_corrected
            result['background'] = m
            result['genes'] = list(sampled)
            results[idx] = result
            idx = idx + 1

        
        # Calculate q-values
        qvals = qvalue.estimate(pvals)

        filtered_results = []

        idx = 0

        for term in sample_map:
            qv = qvals[idx]
            res = results[idx]

            pv = res['p-value']
            k = len(res['genes'])
            res['q-value'] = qv
            if pv < cutoff and k >= gene_threshold:
                filtered_results.append(res)
       
            idx = idx + 1

        return {'results':filtered_results, 'total_genes':total_num_genes}


    def __calculateTermFrequency(self, genes_of_interest, gene2terms):
        sample_term_map = {}

        for gene in genes_of_interest:
            terms = []
            if gene in gene2terms:
                terms = gene2terms[gene]
            else:
                continue
            
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


class TermEnrichment(restful.Resource):
    def get(self):
        return {'results':'use POST to get results.'}

    def post(self):
        results = []
        args = parser.parse_args()

        # Gene names are case INSENSITIVE!
        genes = args[GENES]
        cutoff = args[ALPHA]

        gene_list = genes.split()
        genes_upper = []
        for gene in gene_list:
            genes_upper.append(gene.upper())

        print(genes_upper)
        return tester.performTest(genes_upper, cutoff)

api.add_resource(TermEnrichment, '/enrich')

# Testers
tester = HypergeometricTest()


# Main
if __name__ == '__main__':
    app.run(host='0.0.0.0',debug=True)

