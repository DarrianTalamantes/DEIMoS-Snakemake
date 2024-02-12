import deimos
import numpy as np

 

# Reload previously saved data, excluding scanid column

ms1 = deimos.load(

    "2022_05_10_QC_05.h5",

    key="ms1",

    columns=["mz", "drift_time", "retention_time", "intensity"])