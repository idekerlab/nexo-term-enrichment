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
sgdTable = open("sgd2symbols.txt", "r")
out = open("goterm2genes_2.txt", "w")
out2 = open("gene2goterms_2.txt", "w")

sgd2symbol = {}
sgd2orf = {}
for line in sgdTable:
    l = line.rstrip()
    parts = l.split()
    sgd2symbol[parts[0]] = parts[1]
    if len(parts) == 3:
        orf = parts[2].split('|')[0]
        print(orf)
        sgd2orf[parts[0]] = orf

in_files = glob.glob('*.gene_term')
f = fileinput.input(in_files)

for line in f:
    process_line(line)

f.close()

idmap = {}
for id in id2terms:
    idmap[id] = id2terms[id]
    idmap[sgd2symbol[id]] = id2terms[id]
    if id in sgd2orf:
        idmap[sgd2orf[id]] = id2terms[id]
        print('FOUND######### ' + id)
    else:
        print('Not found: ' +  id)

keys = idmap.keys()
keys.sort()
for id in keys:
    line = id + '\t'
    terms = list(idmap[id])
    terms.sort()
    for term in terms:
        line = line + term + '|'
    
    out2.write(line[0:len(line)-1] + '\n')
    

for term in goterm2symbols:
    sgdIDs = list(goterm2symbols[term])

    line = term + '\t'
    symbols = []
    for sgd in sgdIDs:
        symbols.append(sgd2symbol[sgd])
    symbols.sort()

    for symbol in symbols:
        line = line + symbol + '|'
    out.write(line[0:len(line)-1] + '\n')

out.close()
out2.close()
