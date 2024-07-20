import plotSWOTmap
import numpy as np

plotSWOTmap.plotCalvertMissionMap(logdir='/Users/cproof/processing/deployments/dfo-k999/dfo-k999-20230516/logs/', latlim=[47., 53], lonlim=[-140, -124], start=np.datetime64('2023-05-16'))