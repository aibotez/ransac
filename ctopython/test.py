import ctopytest
import threading,time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
lst = np.zeros((40,10))
lst = lst.tolist()
# print(type(lst))
# lst = [[1,2],[2.5]]
# print(lst[0][0])
s=ctopytest.Test()
# ls=s.Add(10,20)
# print(ls,'lss')




pro_data='86787.txt'
df_news = pd.read_table(pro_data,header = None)
R = df_news[0]
R=np.array(R)
R0 = np.c_[R,R]
# print(R0)
ne = df_news[1]
ne = np.array(ne)
ne0 = np.c_[ne,ne]
for i in range(50):
    R0 = np.c_[R0,R]
    ne0 = np.c_[ne0, ne]


with open('R0.txt', 'a') as f4:
    np.savetxt(f4, R0, delimiter='\t', newline='\n')
f4.close()
with open('ne0.txt', 'a') as f4:
    np.savetxt(f4, ne0, delimiter='\t', newline='\n')
f4.close()


R0 = R0[:,0].tolist()
ne0 = ne0[:,0].tolist()
# print(ne0)
# print(R0)
def run1():
    while True:
        print('cur:', s.GetCurProcess())
def run():
    a = s.ExcuteRan(R0.tolist(), ne0.tolist())
t1 = time.time()
a = s.ExcuteRan(R0, ne0)
npa = np.array(a)
print(np.shape(npa))
# ths = threading.Thread(target=run)
# ths.setDaemon(True)
# ths.start()
# ths1 = threading.Thread(target=run1)
# ths1.setDaemon(True)
# ths1.start()
print(a)
print(time.time()-t1)



