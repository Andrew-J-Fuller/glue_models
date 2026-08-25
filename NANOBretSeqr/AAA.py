import os
os.environ.pop("HTTP_PROXY", None)
os.environ.pop("HTTPS_PROXY", None)
import csv
import time
import requests
from pathlib import Path
import re
import statistics as stats
#  ====================================================================================
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{}"
TLS = r"YOUR PATH HERE" # path to TLS to access uniprot on restricted networks
datafoldr = Path(r'YOUR PATH HERE') # path to folder where all of your exported CSV files are from PRISM
outputfoldr = Path(r'YOUR PATH HERE') # path to folder where you'd like your new formated data
# maybe text wrap?
# when the tag order is "tag"-E3 Ligase
halocrbn = "MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISGEPTTEDLYFQSGSSGMAGEGDQQDAAHNMGNHLPLLPAESEEEDEMEVEDQDSKEAKKPNIINFDTSLPTSHTYLGADMEEFHGRTLHDDDSCQVIPVLPQVMMILIPGQTLPLQLFHPQEVSMVRNLIQKDRTFAVLAYSNVQEREAQFGTTAEIYAYREEQDFGIEIVKVKAIGRQRFKVLELRTQSDGIQQAKVQILPECVLPSTMSAVQLESLNKCQIFPSKPVSREDQCSYKWWQKYQKRKFHCANLTSWPRWLYSLYDAETLMDRIKKQLREWDENLKDDSLPSNPIDFSYRVAACLPIDDVLRIQLLKIGSAIQRLRCELDIMNKCTSLCCKQCQETEITTKNEIFSLSLCGPMAAYVNPHGYVHETLTVYKACNLNLIGRPSTEHSWFPGYAWTVAQCKICASHIGWKFTATKKDMSPQKFWGLTRSALLPTIPDTEDEISPDKVILCLG"
# when the tag is E3 Ligase-"tag"
crbnhalo = "MAGEGDQQDAAHNMGNHLPLLPAESEEEDEMEVEDQDSKEAKKPNIINFDTSLPTSHTYLGADMEEFHGRTLHDDDSCQVIPVLPQVMMILIPGQTLPLQLFHPQEVSMVRNLIQKDRTFAVLAYSNVQEREAQFGTTAEIYAYREEQDFGIEIVKVKAIGRQRFKVLELRTQSDGIQQAKVQILPECVLPSTMSAVQLESLNKCQIFPSKPVSREDQCSYKWWQKYQKRKFHCANLTSWPRWLYSLYDAETLMDRIKKQLREWDENLKDDSLPSNPIDFSYRVAACLPIDDVLRIQLLKIGSAIQRLRCELDIMNKCTSLCCKQCQETEITTKNEIFSLSLCGPMAAYVNPHGYVHETLTVYKACNLNLIGRPSTEHSWFPGYAWTVAQCKICASHIGWKFTATKKDMSPQKFWGLTRSALLPTIPDTEDEISPDKVILCLGSSGEDLYFQSDNDGSEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG"
# below is the nanoluc sequence
nanoluc = "MWLVSLAIVTACAGAMAVYPYDVPDYAGYPYDVPDYAGSYPYDVPDYAGSGVFTLEDFVGDWRQTAGYNLDQVLEQGGVSSLFQNLGVSVTPIQRIVLSGENGLKIDIHVII PYEGLSGDQMGQIEKIFKVVYPVDDHHFKVILHYGTLVIDGVTPNMIDYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLINPDGSLLFRVTINGVTGWRLCERILA"


# uses the protein's common name and the uniprot API to find the protein's AA sequence
# note the protein's common name MUST be the file name (eg FOXA1.csv)
def search_uniprot(protein:str) -> str | None:
    # query parameters to find uniprotID
    params = {
        "query": protein,
        "fields": "accession",
        "size": 1,
    }
    headers = {
        "Accept":"application/json"
    }
        # since some firewalls are pretty picky, we give it a few attempts to reach uniprot
    for attempt in range(1,4):
        try:
            response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=60, verify=TLS)
            results = response.json().get("results",[]) 
            response.raise_for_status()
            if not results:
                return None
            uniprotID =  results[0].get('primaryAccession') # this the uniprot id for the protein
            AA_sequence = requests.get(UNIPROT_ENTRY_URL.format(uniprotID), headers=headers, timeout=60, verify=TLS).json().get('sequence',[]).get('value',[])
            print("Sequence Aquired!")
            return AA_sequence
        
        except (requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout) as e:
            print(f"[attempt {attempt}] error: {e}")
            if attempt == 3:
                raise
            time.sleep(3.0)

# uses search_uniprot, file names, and the above file paths/tag sequences to create a csv shape (#ofproteins,3) where the columns are:
# protein name, sequnece 1 (tag+E3 ligase), sequence 2 (tag + POI), and the mean NANOBret signal resulting
def amino_acid_aquirer(folder:str) -> csv:
    first_sequences = []
    second_sequneces = []
    signal_means = []
    protein_names = []

    for file in folder.iterdir(): # for each file in your datafoldr
        protein = file.name.strip('.csv') 
        AA_sequence = search_uniprot(protein) #aquires the sequence
        time.sleep(0.5)
        protein_names.append(protein)

        with file.open(newline='') as sequences, file.open(newline='') as signals:
            seqs = csv.reader(sequences) 
            sigs = csv.reader(signals)

            # below determines the tag order, and then merges the sequences as would be seen _in vivo_

            for seq in seqs:
              bit = (re.split(r"[-+]",seq[0]))
              first = [] # the E3 ligase and the halo tag combo
              second = [] # the nanoluc tag and the POI combo


              if AA_sequence == None: # if we can't find the uniprotID
                 continue
              if bit == ['ï»¿']: # for some reason empty cells return as /this/
                 continue
              if "Neg" in bit[:2]: # remove the negative controls
                    continue
              if bit[:2] == ["Halo", "CRBN"]:
                  first.append(halocrbn)
              if bit[:2] == ["CRBN", "Halo"]:
                  first.append(crbnhalo)
              if bit[2:] == ["NL", protein]:
                   second = nanoluc + AA_sequence
              if bit[2:] == [protein, "NL"]:
                   second = AA_sequence + nanoluc
              first_sequences.append(first[0])
              second_sequneces.append(second)

              # now we only take the values from the non-negative controls and header (skip 3 rows)

            next(sigs,None)
            next(sigs,None)
            next(sigs,None)
            for sig in sigs: 
                if AA_sequence is None: # if we can't find the uniprotID
                    continue
                values = sig[1:]
                mean = stats.mean(float(value) for value in values)
                signal_means.append(mean)

    # now we write the first sequence, second sequence, and the mean of their nanobret signals to a csv

    with outputfoldr.open('w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for protname,sequence1,sequence2, signal in zip(protein_names, first_sequences, second_sequneces, signal_means):
            writer.writerow([protname,sequence1,sequence2,signal])

amino_acid_aquirer(datafoldr)
