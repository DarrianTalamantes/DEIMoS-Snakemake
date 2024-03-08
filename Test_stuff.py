import deimos
import numpy as np
import matplotlib.pyplot as plt
import pickle 


ms1 = deimos.load(
    "output/samples/parsed/2022_05_10_NR46543_02.h5",
    key="ms1",
    columns=["mz", "drift_time", "retention_time", "intensity"]
)


# Add column names
column_names = ["mz", "drift_time", "retention_time", "intensity"]
header = f"{column_names}\n"


# save data to numpy data frame
np.savetxt('delete_me.csv', ms1, delimiter=',', header=header)
