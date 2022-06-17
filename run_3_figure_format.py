import pickle
import os
import numpy as np
import matplotlib.pyplot as plt

with open('saved_data' + os.sep + 'deserialised.pickle', 'rb') as handle:
    d = pickle.load(handle)

s = d['deserialised']
v = d['vocabulary']
rows = d['rows']
columns = d['columns']
min_pitch = d['min_pitch']
max_pitch = d['max_pitch']

print(s.shape)
print(v.shape)

idx = np.random.randint( s.shape[0] )