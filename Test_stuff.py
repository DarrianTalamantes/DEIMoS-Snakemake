import deimos
import numpy as np
import matplotlib.pyplot as plt
import pickle 


ms1 = deimos.load(
    "output/samples/parsed/2022_05_10_NR46543_02.h5",
    key="ms1",
    columns=["mz", "drift_time", "retention_time", "intensity"]
)
# Convert dictionary values to a NumPy array
data_array = np.array(list(ms1.values()))

# Save the NumPy array to a csv file
np.savetxt("delete_me.csv", data_array, delimiter=',')
