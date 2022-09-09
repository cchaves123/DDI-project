from collections import defaultdict
from collections import Counter
import networkx as nx
import nltk
from nltk.corpus import stopwords
import numpy as np
import pandas as pd
import re
import sys, os


current_dir = os. getcwd()
sys.path.append(current_dir)
from nltk.stem.wordnet import WordNetLemmatizer
from parse_functions import *

nltk.download('omw-1.4')
nltk.download('stopwords')
NAME_FILE = "/Users/cdonnat/Downloads/full database.xml"



file1 = open(NAME_FILE, "r")
str1=file1.read()
file1.close()
list1=str1.split('drug type')

##### Process the data
x=list1[1]
dictionary = {}
for i in range(1,len(list1)):
    x = list1[i]
    res = find_eos(x)
    index_name = res['drugbank_id']
    dictionary[index_name] = res


#trying to find all the most likely occurance word
#Here we didn't run all the 'name' and 'class' key since it is meaningfulless to run them
des=''
ind=''
pha=''
mec=''
tox=''
met=''
abso=''
rou=''
vol=''

for index_name in dictionary.keys():
    newdes = dictionary[index_name]['description']
    des += newdes
    newind = dictionary[index_name]['indication']
    ind += newind
    newpha= dictionary[index_name]['pharmacodynamics']
    pha +=  newpha
    newmec = dictionary[index_name]['mechanism_of_action']
    mec += newmec
    newtox = dictionary[index_name]['toxicity']
    tox += newtox
    newmet = dictionary[index_name]['metabolism']
    met += str(newmet)
    newabso = dictionary[index_name]['absorption']
    abso += newabso
    #newhal = dictionary[index_name]['half_life']
    #hal=hal+""+newhal
    newrou = dictionary[index_name]['route_of_elimination']
    rou += newrou

des= re.split('\W+', des.lower())
fdistdes = Counter(des)
ind= re.split('\W+', ind.lower())
fdistind = Counter(ind)
pha= re.split('\W+', pha.lower())
fdistpha = Counter(pha)
mec=re.split('\W+', mec.lower())
fdistmec = Counter(mec)
tox=re.split('\W+', tox.lower())
fdisttox = Counter(tox)
abso=re.split('\W+', abso.lower())
fdistabso = Counter(abso)

total = fdistdes + fdistind + fdistpha + fdisttox + fdistabso

Lem = WordNetLemmatizer()
wordlistdes=list(total.keys())
invalid_words=list(stop_words)
pattern0=r'^[0-9]+[0-9]$'
regex=re.compile(pattern0)
dict_words =  ['organic','ring', 'acid', 'benzene', 'carbon','alpha', 'fda','aromatic','activity','trial',
               'derivatives',
               'amino', 'moiety', 'allergenic', 'pollen', 'nitrogen', ' cancer',
               'approved', 'plant','drug', 'phenyl', 'bond',
                'pyrimidine', 'receptor', 'protein', 'aromatic', 'allergenic', 'cells','receptor', 'human']
for word in wordlistdes:
    if fdistdes[word]>800 or fdistdes[word]<=1 or h!=None or len(word)<3:
         invalid_words.append(word)

    else:
        w  = Lem.lemmatize(word)
        if w in dict_words:
            pass
        else:
            dict_words += [w]
