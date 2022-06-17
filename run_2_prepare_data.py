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

arr2idx = { str(u):i for i,u in enumerate(v.T) }
idx2arr = { i:u for i,u in enumerate(v.T) }

# add start and end index
arr2idx[str(-1*np.ones( v.shape[0] ))] = -1
idx2arr[-1] = -1*np.ones( v.shape[0] )

print( len(arr2idx) )
print( len(idx2arr) )

in_rows = []
out_rows = []

for i in range(0, s.shape[1]-columns-1, columns):
    inps = np.c_[ -1*np.ones( v.shape[0] ) , s[:, i:i+columns]]
    targ = np.c_[ s[:, i+1:i+columns+1] , -1*np.ones( v.shape[0] ) ]
    in_rows.append( [arr2idx[ str(r) ] for r in inps.T] )
    out_rows.append( [arr2idx[ str(r) ] for r in targ.T] )

print('in_rows:', len(in_rows))
print('out_rows:', len(out_rows))