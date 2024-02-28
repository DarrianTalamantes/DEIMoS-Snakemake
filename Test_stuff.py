import deimos
import numpy as np
import matplotlib.pyplot as plt

ms1 = deimos.load(
    "2022_05_10_NR46543_02.h5",
    key="ms1",
    columns=["mz", "drift_time", "retention_time", "intensity"]
)
print(ms1)
