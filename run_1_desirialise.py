# %%
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt

# %%
with open('saved_data' + os.sep + 'data_tower.pickle', 'rb') as handle:
    d = pickle.load(handle)

# %%
# plot examples
s = d['serialised_segments']
rows = d['rows']
columns = d['columns']
min_pitch = d['min_pitch']
max_pitch = d['max_pitch']

# %%
# deserialize
deserialised = np.zeros( (rows, columns*s.shape[0]) )
tmp_col = 0
for i in range(s.shape[0]):
    deserialised[: ,tmp_col*columns:(tmp_col+1)*columns] = np.reshape( s[i,:], (rows,columns) )
    tmp_col += 1

# %%
plt.imshow(deserialised[:,:200])

# %%
ur = np.unique( deserialised, axis=1 )

# %%
plt.imshow(ur[:,1000:1200])

# %%
with open('saved_data' + os.sep + 'deserialised.pickle', 'wb') as handle:
    d = {}
    d['deserialised'] = deserialised
    d['vocabulary'] = ur
    d['rows'] = rows
    d['columns'] = columns
    d['min_pitch'] = min_pitch
    d['max_pitch'] = max_pitch
    pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)


