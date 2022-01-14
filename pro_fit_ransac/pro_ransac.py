import time,os,random
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from math import *
import pandas as pd
import math,threading

from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor


pro_data='86787.txt'
df_news = pd.read_table(pro_data,header = None)
R = df_news[0]
R=np.array(R)

ne = df_news[1]
ne = np.array(ne)

class ransac():

    def __init__(self,x,y):
        self.x=x
        self.y=y
        self.inside_dist=0.02
        self.C=0
        self.thred_iternum=0
        # self.lock = threading.Lock()
    def reC(self,distance):
        cout_dist = sum(distance)
        insidenums = len([i for i in distance if i <= self.inside_dist])
        # print('27',cout_dist,insidenums)
        C = insidenums/cout_dist
        return C

    def distance_p(self,fit_x,fit_y,point_x,point_y):
        distance=[]
        # ls = [i for i in fit_x if i<= point_x]
        # print(fit_x[fit_x<=point_x])
        # print(fit_x.size)

        idx=np.argwhere(fit_x<=point_x)
        # print(idx)
        if idx.size == 0:
            fitx = fit_x[-1]
            fity = fit_y[-1]
        else:
            idx=idx[0][0]
            fitx = fit_x[idx]
            fity = fit_y[idx]
        # print(point_x,fitx)

        min_distance = ((fitx-point_x)**2+(fity-point_y)**2)**0.5



        # for i in range(len(fit_x)):
        #     distance.append(((fit_x[i]-point_x)**2+(fit_y[i]-point_y)**2)**0.5)
        # min_distance = min(distance)
        return min_distance
    def random_point(self,num):
        # self.lock.acquire()
        if num > len(self.x):
            num = len(self.x)
        data_inx = list(range(1,len(self.x)))
        idx=random.sample(data_inx,num)
        data_x = self.x[idx]
        data_y = self.y[idx]
        # self.lock.release()

        return [data_x,data_y]

    def liner_fit(self,x, A1, B1):
        return A1 * x + B1
    def mtanh_2(self,x,A0,B,alpha,x_sys,w,beta):
        z = (x_sys - x)/ w
        # for i in range(len(z)):
        #     if abs(z[i]) >500:
        #         z[i] =500
        mtanh_fun = ((1 + alpha * z)* np.exp(z) - (1 + beta * z)* np.exp(-z))/ (np.exp(z) + np.exp(-z))
        y = A0 * mtanh_fun + B
        return y
    def fit_act1(self,x,y,nh_type):
        if nh_type == 'L':
            A1, B1 = optimize.curve_fit(self.liner_fit, x, y)[0]
        elif nh_type=='M':
            ped_pos = 2.305 # pedestal position
            h = 2 # height
            w = 0.02 # width
            slope1 = 0.1 # slope of the core part
            slope2 = -0.07 # slope of edge part
            a0 = [h / 2,h / 2,slope1,ped_pos,w,slope2]
            A,B = optimize.curve_fit(self.mtanh_2,x,y,a0,maxfev=1000)
            fit_x = np.arange(self.x[0],self.x[-1], -0.001)  # 30和75要对应x0的两个端点，0.01为步长
            fit_y = self.mtanh_2(fit_x, A[0], A[1],A[2],A[3],A[4],A[5])
        fit_data = [fit_x, fit_y]
        return fit_data
    def fit_act(self,x,y,nh_type):
        if nh_type == 'L':
            A1, B1 = optimize.curve_fit(self.liner_fit, x, y)[0]
        elif nh_type=='M':
            ped_pos = 2.305 # pedestal position
            h = 2 # height
            w = 0.02 # width
            slope1 = 0.1 # slope of the core part
            slope2 = -0.07 # slope of edge part
            a0 = [h / 2,h / 2,slope1,ped_pos,w,slope2]
            A,B = optimize.curve_fit(self.mtanh_2,x,y,a0,maxfev=200)
            fit_x = np.arange(self.x[0],self.x[-1], -0.001)  # 30和75要对应x0的两个端点，0.01为步长
            fit_y = self.mtanh_2(fit_x, A[0], A[1],A[2],A[3],A[4],A[5])
        distance_list=[]
        for i in range(len(self.x)):
            jl = self.distance_p(fit_x,fit_y,self.x[i],self.y[i])
            distance_list.append(jl)
        # plt.subplot(211)
        # plt.scatter(self.x,self.y)
        # plt.plot(fit_x, fit_y)
        # plt.subplot(212)
        # plt.scatter(self.x,distance_list)
        # plt.show()
        fit_data = [fit_x,fit_y]
        return [distance_list,fit_data]

    def thread_act(self):
        num=self.thred_iternum
        C=[]
        data_fit=[]
        for i in range(num):
            # print(i)
            rand_data = self.random_point(6)
            rand_x = rand_data[0]
            rand_y = rand_data[1]
            try:
                fit_re=self.fit_act(rand_x,rand_y,'M')
                dist = fit_re[0]
                rec = self.reC(dist)
                if not math.isnan(rec):
                    C.append(self.reC(dist))
                    data_fit.append(fit_re[1])
            except:
                pass
        C_max = max(C)
        idx = C.index(C_max)
        data = data_fit[idx]
        return [C_max,data]

    def iter_act(self,num):
        t1 = time.time()
        thrednum=5
        self.thred_iternum = int(num/thrednum)
        executor = ProcessPoolExecutor(max_workers=thrednum)
        futures = []
        for i in range(0,thrednum):
            future = executor.submit(self.thread_act)
            futures.append(future)
        executor.shutdown(True)
        t2 = time.time()
        print('use_time: ',t2-t1)
        C=[]
        data = []
        for future in futures:
            result = future.result()
            C.append(result[0])
            data.append(result[1])
        C_max = max(C)
        print('Cmax', C_max)
        idx = C.index(C_max)
        data0 = data[idx]
        st1=time.time()
        for ri in range(1):
            fit_re = self.fit_act(self.x, self.y, 'M')
        print(time.time()-st1)
        fit_x = fit_re[1][0]
        fit_y = fit_re[1][1]

        plt.scatter(self.x, self.y)
        plt.plot(data0[0], data0[1])
        plt.plot(fit_x, fit_y, 'blue')
        plt.show()





        # C=[]
        # data_fit = []
        # for i in range(num):
        #     # print(i)
        #     rand_data = self.random_point(6)
        #     rand_x = rand_data[0]
        #     rand_y = rand_data[1]
        #     try:
        #         fit_re = self.fit_act(rand_x, rand_y, 'M')
        #         dist = fit_re[0]
        #         rec = self.reC(dist)
        #         if not math.isnan(rec):
        #             C.append(self.reC(dist))
        #             data_fit.append(fit_re[1])
        #     except:
        #         pass
        # t2 = time.time()
        # print('use_time: ', t2 - t1)
        # C_max = max(C)
        # print('Cmax',C_max)
        # idx = C.index(C_max)
        # data = data_fit[idx]
        # fit_re = self.fit_act(self.x, self.y, 'M')
        # fit_x=fit_re[1][0]
        # fit_y = fit_re[1][1]
        # plt.scatter(self.x, self.y)
        # plt.plot(fit_x, fit_y,'blue')
        # plt.plot(data[0], data[1])
        # plt.show()

        return C_max
if __name__ == '__main__':
    a=ransac(R,ne)
    a.iter_act(3000)
# rand_data=a.random_point(6)
# rand_x = rand_data[0]
# rand_y = rand_data[1]
# dis=a.fit_act(rand_x,rand_y,'M')
















