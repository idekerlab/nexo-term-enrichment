#!/opt/local/bin/python
import json

# Input data files
nodeTable = open("nexo_original.json", "r")
out = open("term2genes.txt", "w")

# Load entire data
nexo3 = json.load(nodeTable)
nodeTable.close()
nexo3Nodes = nexo3["elements"]["nodes"]

gene2term = {}
symbol2sgd = {}

def addTerms(genes):
  for gene in genes:
    if gene in gene2term:
      termList = gene2term[gene]
      termList.append(termID)
      gene2term[gene] = termList
    else:
      termList = []
      termList.append(termID)
      gene2term[gene] = termList


for node in nexo3Nodes:
  termID = node["data"]["NeXO Term ID / SGD Gene ID"]
  geneNames = node["data"]["Assigned Genes"]
  orfNames = node["data"]["Assigned Orfs"]

  if termID.startswith("S"):
    geneNames = geneNames.replace("'", "")
    geneNames = geneNames.replace("\"", "")
    symbol2sgd[geneNames.split(",")[0]] = termID
  else:
    geneNames = geneNames.replace("]", "")
    geneNames = geneNames.replace("[", "")
    geneNames = geneNames.replace("'", "")
    geneNames = geneNames.replace("\"", "")
    geneNames = geneNames.replace(" ", "")
    orfNames = orfNames.replace("'","")
    orfNames = orfNames.replace("]","")
    orfNames = orfNames.replace("[","")
    orfNames = orfNames.replace(" ","")
    genes = geneNames.split(",")
    orfs = orfNames.split(",")
    addTerms(genes)
    addTerms(orfs)
    out.write(termID + "\t" + geneNames + "\n")

for key in symbol2sgd:
  if key in gene2term:
    terms = gene2term[key]
    gene2term[symbol2sgd[key]] = terms
  else:
    print("Error! ===> " + key)

for key in gene2term:
  terms = gene2term[key]
  terms.sort()
  termStr = ""
  for  term in terms:
    termStr = termStr + term + ","

  print(key + "\t" + termStr[0:len(termStr)-1])

out.close()

