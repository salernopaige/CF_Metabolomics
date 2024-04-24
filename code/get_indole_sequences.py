from Bio import Entrez

Entrez.email = "salernopaige@gmail.com"
Entrez.api_key = "c8e086936ec0e50f4f721a9c37ae4b546e08"

search_handle = Entrez.esearch(db="nucleotide", term = "tnaA bacteria", retmax=1000)
search_results = Entrez.read(search_handle)
search_handle.close()

ids = search_results["IdList"]
fetch_handle = Entrez.efetch(db = "nucleotide", id=ids, rettype="fasta", retmode="text")
data=fetch_handle.read()
fetch_handle.close()

with open("tnaA_sequences.fasta", "w") as fasta_file:
    fasta_file.write(data)
