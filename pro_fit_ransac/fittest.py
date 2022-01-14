import time,os,random
from numba import cuda # 从numba调用cuda
cuda.select_device(0)
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from math import *
import pandas as pd
from scipy import interpolate
import math,threading
from scipy.optimize import curve_fit
from numba import jit

pro_data='86787.txt'
df_news = pd.read_table(pro_data,header = None)
R = df_news[0]
x=np.array(R)

ne = df_news[1]
y = np.array(ne)
xnew=np.linspace(2.3,2.278,100)
f = interpolate.interp1d(x, y, kind='slinear')
ynew = f(xnew)

xnew = x
ynew = y

def func1(x,a,b):
    return a*np.exp(b*x) - a


@jit(nopython=True, cache=True)
def func(x):
    coeff = np.polyfit(x,x, 1)
    return exp(x)

@cuda.jit
def Fittest(N,x,y):
    idx = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x
    if (idx < N):
        # 拟合
        yy = func(x[idx])
        # popt, pcov = curve_fit(func, x, y)


#生成x,y数据
x = np.linspace(0.6,1.7,23)
y = func1(x,0.5,2.0)
y = y + 0.1 * np.random.randn(len(x))

XG = cuda.to_device(xnew)
YG = cuda.to_device(ynew)
Fittest[32, 100](10,XG,YG)
cuda.synchronize()




#
# a=np.polyfit(xnew,ynew,3)#用2次多项式拟合x，y数组
# b=np.poly1d(a)#拟合完之后用这个函数来生成多项式对象
# c=b(xnew)#生成多项式对象之后，就是获取x在这个多项式处的值
# plt.scatter(xnew,ynew,marker='o',label='original datas')#对原始数据画散点图
# plt.plot(xnew,c,ls='--',c='red',label='fitting with second-degree polynomial')#对拟合之后的数据，也就是x，c数组画图
# plt.legend()
# plt.show()