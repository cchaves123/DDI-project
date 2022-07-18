from collections import defaultdict
import networkx as nx
import numpy as np
import pandas as pd
import re
import sys, os


current_dir = os. getcwd()
sys.path.append(current_dir)
from parse_functions import *

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
