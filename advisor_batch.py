import os
import json
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors

# Number of results. Can be changed
length = "15"

affinity = input("Enter affinity (strong, standard, or weak): ")

read_smiles = open("smiles.txt", "r")
query_results = open("results.txt", "w")
term_size = os.get_terminal_size()
query_results.write('=' * term_size.columns + "\n")

for i in read_smiles:

    # Get request Arthor Similarity Search Data
    url = "https://arthor.docking.org/dt/aggpage/search?query=" + i.replace('\n', '') + "&type=Similarity&start=0&length=" + length + "&draw=0&fmt=json"
    response = requests.get(url)
    result = response.json()

    query_logP = round(Descriptors.MolLogP(Chem.MolFromSmiles(str(result["query"]))), 3)
    message = ""
    similarity = float(result["data"][0][2])

    if (similarity == 1.000 ):
        message = "This molecule has previously been observed to aggregate, as below."
        if (affinity == "strong"):
            message += " At this high affinity, the activity may be genuine. We recommend performing experimental controls."
    elif (similarity < 0.450):
        message = "This molecule is very similar to a molecule that has been previously observed to aggregate."
        if (affinity == "strong"):
            message += " At this high affinity, the activity may be genuine. We recommend performing experimental controls."
    elif (similarity > 0.400 and query_logP > 3.500):
        message = "This molecule is somewhat similar to a molecule that has been previously been observed to aggregate, below, AND, with a high logP, is suspicious. We recommend performing experimental controls."
        if (affinity == "standard" or affinity == "weak"):
            message += " Controls are particularly important to run at micromolar concentrations."
    elif (similarity <= 0.400 and query_logP > 3.500):
        message = "This molecule is not similar to a known aggregator, but with a high logP, is a candidate for aggregation. We recommend performing experimental controls."
    else:
        message = "This molecule does not look like one that has been previously observed to aggregate. Although the risk is low based on past experience, we always recommend performing experimental controls to confirm the lack of artifactual inhibition due to aggregation."

    query_results.write("Query: {}\nLogP: {}\n{}\n{}\n".format(str(result["query"]), query_logP, message,('=' * term_size.columns)))


    for j in result["data"]:
        temp = j[1].split("~")
        long_string = (temp[0].split(" "))
        smiles = long_string[0]
        smile_name = long_string[1]
        smile_logP = round(Descriptors.MolLogP(Chem.MolFromSmiles(str(smiles))), 3)
        query_results.write("{:<5} {:<20} \t\t {:<20}{:<20} {:<10} \n".format(str(j[0]), str(smiles), str(smile_name), str(j[2]), smile_logP))
    query_results.write("\n" + '=' * term_size.columns + "\n")

query_results.close()
read_smiles.close()
