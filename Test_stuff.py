import deimos
import numpy as np
import matplotlib.pyplot as plt
import pickle 


ms1 = deimos.load(
    "2022_05_10_NR46543_02.h5",
    key="ms1",
    columns=["mz", "drift_time", "retention_time", "intensity"]
)
# Need to save this dictionary not print it
print(ms1)
with open('delete_me.pkl', 'wb') as f:
    pickle.dump(ms1, f)
np.save('delete_me.npy', ms1)
