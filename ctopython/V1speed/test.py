import ctopytest
import threading,time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

s=ctopytest.Test()



pro_data='86787.txt'
R  = pd.read_table('R.txt',header = None)
R=np.array(R)
[row,col] = np.shape(R)
R = R.tolist()
ne = pd.read_table('ne.txt',header = None)
ne = np.array(ne).tolist()




t1 = time.time()
a = s.ExcuteRan(R, ne)

print(a)
print(time.time()-t1)



