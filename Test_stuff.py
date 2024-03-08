import deimos
import numpy as np
import matplotlib.pyplot as plt
import pickle 


ms1 = deimos.load(
    "output/samples/parsed/2022_05_10_NR46543_02.h5",
    key="ms1",
    columns=["mz", "drift_time", "retention_time", "intensity"]
)
# save data to numpy data frame
np.savetxt('delete_me.npy', ms1, delimiter=',')
