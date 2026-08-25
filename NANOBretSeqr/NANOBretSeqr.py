''Sources:
'https://pubmed.ncbi.nlm.nih.gov/18247345/' (Bandopadhyay-Mehler hydrophobicity scale | See Table 2)

'''
import os

os.environ.pop("HTTP_PROXY", None)
os.environ.pop("HTTPS_PROXY", None)

import numpy as np
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import requests
from pathlib import Path
import csv
import time
import re
import statistics as stats

#=========================  CREATE TRAINING DATA  ================================================================
#=================================================================================================================

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{}"
TLS = r"YOUR PATH HERE"
datafoldr = Path(r"YOUR PATH HERE")
trainingfoldr = Path(r'YOUR PATH HERE')
halocrbn = "MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISGEPTTEDLYFQSGSSGMAGEGDQQDAAHNMGNHLPLLPAESEEEDEMEVEDQDSKEAKKPNIINFDTSLPTSHTYLGADMEEFHGRTLHDDDSCQVIPVLPQVMMILIPGQTLPLQLFHPQEVSMVRNLIQKDRTFAVLAYSNVQEREAQFGTTAEIYAYREEQDFGIEIVKVKAIGRQRFKVLELRTQSDGIQQAKVQILPECVLPSTMSAVQLESLNKCQIFPSKPVSREDQCSYKWWQKYQKRKFHCANLTSWPRWLYSLYDAETLMDRIKKQLREWDENLKDDSLPSNPIDFSYRVAACLPIDDVLRIQLLKIGSAIQRLRCELDIMNKCTSLCCKQCQETEITTKNEIFSLSLCGPMAAYVNPHGYVHETLTVYKACNLNLIGRPSTEHSWFPGYAWTVAQCKICASHIGWKFTATKKDMSPQKFWGLTRSALLPTIPDTEDEISPDKVILCLG"
crbnhalo = "MAGEGDQQDAAHNMGNHLPLLPAESEEEDEMEVEDQDSKEAKKPNIINFDTSLPTSHTYLGADMEEFHGRTLHDDDSCQVIPVLPQVMMILIPGQTLPLQLFHPQEVSMVRNLIQKDRTFAVLAYSNVQEREAQFGTTAEIYAYREEQDFGIEIVKVKAIGRQRFKVLELRTQSDGIQQAKVQILPECVLPSTMSAVQLESLNKCQIFPSKPVSREDQCSYKWWQKYQKRKFHCANLTSWPRWLYSLYDAETLMDRIKKQLREWDENLKDDSLPSNPIDFSYRVAACLPIDDVLRIQLLKIGSAIQRLRCELDIMNKCTSLCCKQCQETEITTKNEIFSLSLCGPMAAYVNPHGYVHETLTVYKACNLNLIGRPSTEHSWFPGYAWTVAQCKICASHIGWKFTATKKDMSPQKFWGLTRSALLPTIPDTEDEISPDKVILCLGSSGEDLYFQSDNDGSEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG"
nanoluc = "MWLVSLAIVTACAGAMAVYPYDVPDYAGYPYDVPDYAGSYPYDVPDYAGSGVFTLEDFVGDWRQTAGYNLDQVLEQGGVSSLFQNLGVSVTPIQRIVLSGENGLKIDIHVII PYEGLSGDQMGQIEKIFKVVYPVDDHHFKVILHYGTLVIDGVTPNMIDYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLINPDGSLLFRVTINGVTGWRLCERILA"

# uses protein common name (acquired from file name) and queries uniprot and returns the ammino acid sequence of the protein
def search_uniprot(protein:str) -> str | None:
    params = {
        "query": protein,
        "fields": "accession",
        "size": 1,
    }
    headers = {
        "Accept":"application/json"
    }
    # tricky firewalls call for persistent requests
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

# using a folder of csvs (eg "FOXA1.csv") returns
# 1. a 2d np.array where each row is a PPI and each column is each protein's AA sequence including the tag
# 2. a 1d np.array where each column is the mean NANOBret signal ratio (averaged from 8 replicates) 
def seqr(folder:str) -> np.array:
    sequencerows = []
    signalrows = []

    for file in folder.iterdir(): 
        protein = file.name.strip('.csv') 
        sequence = search_uniprot(protein) 
        time.sleep(0.5)

            # aquire the sequences
        with file.open(newline='') as csvfile: 
            reader = csv.reader(csvfile)
            for row in reader: 
              bit = (re.split(r"[-+]",row[0]))
              front = [] #the E3 ligase and it's HALO Tag combo
              back = [] #the NanoLuc tag and the POI combo
              if sequence == None: 
                 continue
              if bit == ['ï»¿']:
                 continue
              if "Neg" in bit[:2]:
                    continue
              if bit[:2] == ["Halo", "CRBN"]:
                  front.append(halocrbn)
              if bit[:2] == ["CRBN", "Halo"]:
                  front.append(crbnhalo)
              if bit[2:] == ["NL", protein]:
                   back = nanoluc + sequence
              if bit[2:] == [protein, "NL"]:
                   back = sequence + nanoluc
              rows = np.array([front[0],back])
              sequencerows.append(rows)


        # aquire the signal ratios
        with file.open(newline='') as csvf1:
            reader = csv.reader(csvf1)
            next(reader,None)
            next(reader,None)
            next(reader,None)
            for row in reader:
                if sequence is None:
                  continue
                val = row[1:]
                mean = stats.mean(float(x) for x in val)
                signalrows.append(mean)

    return np.vstack(sequencerows), np.array(signalrows)


#=========================  DEFINE AMINO ACID ENCODING  ================================================================
#=======================================================================================================================

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
PAD_TOKEN = " "
alphabet = PAD_TOKEN + AMINO_ACIDS
char_to_idx = {ch: i for i, ch in enumerate(alphabet)}
max_seq_len = 600 # max length for each protein sequence possible
n_tokens = len(alphabet) #the number of unique codes in the amino acid alphabet (20 standard plus one pad token)

# in an attempt to featurize each amino acid the following properties will be added to each
AA_prop = {
    # the order of property is: 
    # (hydrophobicity (Bandyopadhyay-Mehler hydrophobicity scale), charge(literal), polar(0/1))
    'A': (0.33,0,0), 'C' : (1.15,0,1), 'D' : (-0.22,-1,1), 'E' : (-0.24,-1,1), 
    'F' : (0.85,0,0), 'G' : (0.01,0.0,0), 'H': (0.25,0,1), 'I' : (0.97,0,0), 
    'K' : (-0.40,1,1), 'L' : (0.87,0,0), 'M' : (0.54, 0, 0), 'N' : (-0.07, 0, 1), 
    'P' : (0.32, 0, 0), 'Q' : (-0.05, 0, 1), 'R' : (-0.01, 1, 1), 'S' : (0.05, 0, 1),
    'T' : (0.21, 0, 1), 'V' : (0.83, 0, 0), 'W' : (0.67, 0, 0), 'Y' : (0.60, 0, 1), ' ': (0,0,0)
}
n_props = 3

#=====================  ENCODE AMINO ACID SEQUENCES  ===================================================================
#=======================================================================================================================

# encode the AA sequence using one_hot arrays
def encodeseq(seq: str, max_len: int = max_seq_len) -> np.ndarray:
    # cleans up the passed sequence
    seq = seq.upper().strip()
    seq = ''.join(ch if ch in AMINO_ACIDS else PAD_TOKEN for ch in seq)
    if len(seq) > max_len:
        seq = seq[:max_len]
    else:
        seq = seq + PAD_TOKEN * (max_len-len(seq))

    # a simple one_hot encoding for each amino acid in the sequence
    rows = []
    for ch in seq:
        one_hot = np.zeros(n_tokens, dtype=float)
        one_hot[char_to_idx[ch]] = 1.0 
        props = np.array(AA_prop[ch], dtype = float)
        rows.append(np.concatenate([one_hot, props])) 

    return np.stack(rows, axis=0) # shape: (max_len, n_tokens + n_props)
   
# reduces the amount of featurization to reduce overfit of the model given smaller sample sizes
def agg(array):
        return np.concatenate([array.mean(axis=0), array.sum(axis=0)])

# combines two amino strings into one np.ndarray for the model (with some extra interactions)
def encodepair(seq1: str, seq2: str) -> np.ndarray:

    m1 = encodeseq(seq1)
    m2 = encodeseq(seq2)
    interaction = m1 * m2 
    difference = m1 - m2 
    result = np.concatenate([agg(m1), agg(m2), agg(interaction), agg(difference)]) 
    return result

# call's seqr to find the AA sequence and signal ratios for each data file in "datafoldr"
def maketrainingdata():
    x, y = seqr(datafoldr)
    seqlist = []
    for seq1, seq2 in x:
        seqlist.append(encodepair(seq1,seq2))
    return np.vstack(seqlist), y

#========================  TRAINING THE MODEL  =========================================================================
#=======================================================================================================================

# takes the traing data from above and shoves it into the model
# prints some pretty things and takes the score from the model too
def train_model():
    print("====  Training Model  ====")

    X, y = maketrainingdata()
    xtrain, xtest, ytrain, ytest = train_test_split(X,y,test_size=20)

    model = MLPRegressor(
        hidden_layer_sizes=(1024,512,256,128),
        activation='relu',
        max_iter=1000000,
        random_state=0,
        alpha = 0.01,
        learning_rate_init = 0.000001,
        early_stopping=False,
    )

    model.fit(xtrain,ytrain)
    ypred = model.predict(xtest)
    mse = mean_squared_error(ytest, ypred)
    
    r2 = model.score(xtrain, ytrain)

    print(f"Test MSE: {mse}")
    print(f"Training R^2 score: {r2}")
    return model

#========================  PREDICTING NANOBRET SIGNAL FROM AMINO ACID PAIR  ============================================
#=======================================================================================================================

# asks the model to predict the NANOBret signal of two amino acid sequences
def predict_signal(model: MLPRegressor, seq1: str, seq2: str):
    features = encodepair(seq1, seq2)
    features = features.reshape(1,-1)
    predicted_signal = model.predict(features)[0]

    print(f"\nInput sequences:")
    print(f"  seq1 = {seq1}")
    print(f"  seq2 = {seq2}")
    print(f"Predicted NanoBRET signal: {predicted_signal:.3f}")




if __name__ == "__main__":
    # pip install scikit-learn
    # python Beginner nanobret.py
    model = train_model()
    
    new_seq1 = "MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISGEPTTEDLYFQSGSSGMAGEGDQQDAAHNMGNHLPLLPAESEEEDEMEVEDQDSKEAKKPNIINFDTSLPTSHTYLGADMEEFHGRTLHDDDSCQVIPVLPQVMMILIPGQTLPLQLFHPQEVSMVRNLIQKDRTFAVLAYSNVQEREAQFGTTAEIYAYREEQDFGIEIVKVKAIGRQRFKVLELRTQSDGIQQAKVQILPECVLPSTMSAVQLESLNKCQIFPSKPVSREDQCSYKWWQKYQKRKFHCANLTSWPRWLYSLYDAETLMDRIKKQLREWDENLKDDSLPSNPIDFSYRVAACLPIDDVLRIQLLKIGSAIQRLRCELDIMNKCTSLCCKQCQETEITTKNEIFSLSLCGPMAAYVNPHGYVHETLTVYKACNLNLIGRPSTEHSWFPGYAWTVAQCKICASHIGWKFTATKKDMSPQKFWGLTRSALLPTIPDTEDEISPDKVILCLG" 
    new_seq2 = "MWLVSLAIVTACAGAMAVYPYDVPDYAGYPYDVPDYAGSYPYDVPDYAGSGVFTLEDFVGDWRQTAGYNLDQVLEQGGVSSLFQNLGVSVTPIQRIVLSGENGLKIDIHVII PYEGLSGDQMGQIEKIFKVVYPVDDHHFKVILHYGTLVIDGVTPNMIDYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLINPDGSLLFRVTINGVTGWRLCERILAMAWKSGGASHSELIHNLRKNGIIKTDKVFEVMLATDRSHYAKCNPYMDSPQSIGFQATISAPHMHAYALELLFDQLHEGAKALDVGSGSGILTACFARMVGCTGKVIGIDHIKELVDDSVNNVRKDDPTLLSSGRVQLVVGDGRMGYAEEAPYDAIHVGAAAPVVPQALIDQLKPGGRLILPVGPAGGNQMLEQYDKLQDGSIKMKPLMGVIYVPLTDKEKQWSRWK"

    # above is seq 1 (halocbrn) and seq2 (nanoluc-PIMT_human) a well described target of CBRN and should result in a positive NANOBret signal ratio
    predict_signal(model, new_seq1, new_seq2)
