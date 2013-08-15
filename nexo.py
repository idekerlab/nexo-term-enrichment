# Term Enrichement
import numpy as np
import scipy as sp

import math

import qvalue

from scipy import stats

from flask import Flask
from flask.ext import restful
from flask.ext.restful import Resource, fields, reqparse

import json

# Constants
GENES = 'genes'
ALPHA = 'alpha'
ONTOLOGY_TYPE = 'type'

#root terms
ROOT_NEXO = "joining_root"
ROOT_BP = "GO:0008150"
ROOT_MF = "GO:0003674"
ROOT_CC = "GO:0005575"

# Valid namespace
VALID_ONTOLOGY = ['NEXO', 'MF', 'CC', 'BP']

# p-value cutoff
p_val_threshold = 0.01

# At least this number of genes should be assigned
gene_threshold = 2

app = Flask(__name__)
api = restful.Api(app)

parser = reqparse.RequestParser()
parser.add_argument(GENES, type=str)
parser.add_argument(ALPHA, type=float)
parser.add_argument(ONTOLOGY_TYPE, type=str)


class TermMapping:

    _gene2terms = {}
    _term2genes = {}
    _num_genes = {}

    def __init__(self):
        print "Loading resource files..."
        
        # NeXO resources
        self._gene2terms['NEXO'] = self.createMap("data/gene2terms.txt")
        self._term2genes['NEXO'] = self.createMap("data/term2genes.txt")
       
        # GO resources
        self._gene2terms['GO'] = self.createMap("data/gene2goterms.txt")
        self._term2genes['GO'] = self.createMap("data/goterm2genes.txt")
        
        # Count genes.  Root term has all genes.
        all_genes = self._term2genes['NEXO'][ROOT_NEXO]
        self._num_genes['NEXO'] = len(all_genes)
        
        all_genes = self._term2genes['GO'][ROOT_MF]
        self._num_genes['MF'] = len(all_genes)
        all_genes = self._term2genes['GO'][ROOT_CC]
        self._num_genes['CC'] = len(all_genes)
        all_genes = self._term2genes['GO'][ROOT_BP]
        self._num_genes['BP'] = len(all_genes)

        print("Ready.  NeXO Genes = " + str(self._num_genes['NEXO']))
        print("MF Genes = " + str(self._num_genes['MF']))
        print("BP Genes = " + str(self._num_genes['BP']))
        print("CC Genes = " + str(self._num_genes['CC']))

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


    def get_gene_mapping(self, ontology_type):
        if ontology_type == 'NEXO':
            return self._gene2terms[ontology_type]
        else:
            return self._gene2terms['GO']

    def get_term_mapping(self, ontology_type):
        if ontology_type == 'NEXO':
            return self._term2genes[ontology_type]
        else:
            return self._term2genes['GO']

    def get_gene_count(self, ontology_type):
        return self._num_genes[ontology_type]


class HypergeometricTest:

    _mapper = None

    def __init__(self):
        self._mapper = TermMapping()

    def performTest(self, ontology_type, genes_of_interest, cutoff):

        if cutoff is None:
            cutoff = p_val_threshold

        total_num_genes = self._mapper.get_gene_count(ontology_type)
        term2genes = self._mapper.get_term_mapping(ontology_type)
        gene2terms = self._mapper.get_gene_mapping(ontology_type)

        sample_map = self.__calculateTermFrequency(genes_of_interest, gene2terms)

        n = len(genes_of_interest)

        # Number of tests performed: will be used for correction.
        num_tests = len(sample_map)

        results = [ {} ] * num_tests

        pvals = np.zeros(num_tests)

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

        
        # Correct border values (library does not accept 0 & 1)
        i = 0
        for p in pvals:
            if p >= 1.0 or math.isnan(p):
                pvals[i] = 0.9999999999999999999999
            elif p <= 0:
                pvals[i] = 0.0000000000000000000001
            i += 1
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
        ontology_type = args[ONTOLOGY_TYPE]

        if ontology_type is None:
            ontology_type = 'NEXO'
        elif ontology_type not in VALID_ONTOLOGY:
            return {'results':'Ontology type not supported: ' + ontology_type}

        gene_list = genes.split()
        genes_upper = []
        for gene in gene_list:
            genes_upper.append(gene.upper())

        #print(genes_upper)
        return tester.performTest(ontology_type, genes_upper, cutoff)

api.add_resource(TermEnrichment, '/enrich')

# Testers
tester = HypergeometricTest()


# Main
if __name__ == '__main__':
    app.run(host='0.0.0.0',debug=True)

