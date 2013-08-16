#!/opt/local/bin/python
import sys
import fileinput
import glob

goterm2symbols = {}
id2terms = {}

def process_line(line):
    l = line.rstrip()
    parts = l.split()
    symbols = set([])
    if parts[1] in goterm2symbols:
        symbols = goterm2symbols[parts[1]]
    
    symbols.add(parts[0])
    goterm2symbols[parts[1]] = symbols

    # build gene-to-terms list
    terms = set([])
    if parts[0] in id2terms:
        terms = id2terms[parts[0]]
    terms.add(parts[1])
    id2terms[parts[0]] = terms



# Input data files
interm = open("goterm2genes_2.txt", "r")
ingenes = open("gene2goterms_2.txt", "r")


#########################
# for ID Mapping
##########################
sgdTable = open("sgd2symbols.txt", "r")
sgd2symbol = {}
sgd2orf = {}
symbol2sgd = {}

for line in sgdTable:
    l = line.rstrip()
    parts = l.split()
    sgd2symbol[parts[0]] = parts[1]
    symbol2sgd[parts[1]] = parts[0]
    if len(parts) == 3:
        orf = parts[2].split('|')[0]
        sgd2orf[parts[0]] = orf


term2genes = {}
gene2terms = {}

for line in interm:
    l = line.rstrip()
    parts = l.split()
    genes = parts[1].split('|')
    term2genes[parts[0]] = set(genes)

for line in ingenes:
    l = line.rstrip()
    parts = l.split()
    terms = parts[1].split('|')
    gene2terms[parts[0]] = set(terms)


trees = ["biological_process.info_gain.gene_term", "cellular_component.info_gain.gene_term","molecular_function.info_gain.gene_term"]


def create_files(tree_name, terms, genes):

    geneout = open(tree_file_name + '.genes', "w")
    termout = open(tree_file_name + '.terms', "w")


    assigned = []

    for term in terms:
        assigned_genes = term2genes[term]
        # filter genes not in this tree.
        symbols = set([])
        for symbol in assigned_genes:
            sgd = symbol2sgd[symbol]
            if sgd in genes:
                symbols.add(symbol)
       
        assigned.append(len(symbols))

        term_line = term + '\t'
        for s in symbols:
            term_line = term_line + s + '|'
        termout.write(term_line[0:len(term_line)-1] + '\n')

    print(tree_name + ' MAX = ' + str(max(assigned)))
    print(tree_name + ' ALL = ' + str(len(genes)))

    for sgd in genes:
        assigned_terms = gene2terms[sgd]
        # Filter terms only in the current tree.
        goterms = set([])
        for term in assigned_terms:
            if term in terms:
                goterms.add(term)

        gene_line = sgd + '\t'
        for t in goterms:
            gene_line = gene_line + t + '|'

        sgd_line = gene_line[0:len(gene_line)-1] + '\n'
        geneout.write(sgd_line)
        
        if sgd in sgd2orf:
            orf_line = sgd_line.replace(sgd, sgd2orf[sgd])
            geneout.write(orf_line)
        
        if sgd in sgd2symbol:
            symbol_line = sgd_line.replace(sgd, sgd2symbol[sgd])
            geneout.write(symbol_line)
    
    geneout.close()
    termout.close()



for tree_file_name in trees:
    # Extract unique terms
    termset = set([])
    geneset = set([])
    tree_file = open(tree_file_name, "r")
    for line in tree_file:
        l = line.rstrip()
        parts = l.split()
        term = parts[1]
        gene = parts[0]
        termset.add(term)
        geneset.add(gene)

    tree_file.close()
    
    create_files(tree_file_name, termset, geneset)







