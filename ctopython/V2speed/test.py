import ctopytest
import threading,time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

s=ctopytest.Test()



pro_data='86787.txt'
R  = pd.read_table('R0.txt',header = None)
R=np.array(R)
[row,col] = np.shape(R)
ne = pd.read_table('ne0.txt',header = None)
ne = np.array(ne)




t1 = time.time()

for i in range(col):
    print(i)
    Ri = R[:,i].tolist()
    nei = ne[:,i].tolist()
    a = s.ExcuteRan(Ri, nei,4)

print(a)
print(time.time()-t1)



