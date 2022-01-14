import os,time,requests,json
import threading
import Main
import numpy as np
import pandas as pd
from scipy import io
import matplotlib.pyplot as plt




class CCOper():
    def __init__(self):
        self.MainExc = Main.ExcuMain()
        self.counts = 0

    def SaveDatatoMat(self,pedestal):
        io.savemat('pedestal.mat', {'pedestal': pedestal})

    def JudgeServerRun(self):
        configurationInfo = self.MainExc.configurationInfo
        if configurationInfo[1] == '1':
            res = self.MainExc.SubmitComputerInfos()
            try:
                res = requests.get('http://127.0.0.1:9988/JudgeServerRun').text
            except:
                os.system('start startproject.cmd')
                time.sleep(5)
        else:
            pass

    def showResult(self,pedestal):
        plt.close()
        fig, axes = plt.subplots(3, 1)
        ax1 = axes[0]
        ax2 = axes[1]
        ax3 = axes[2]
        for i in axes:
            i.cla()
        [row, col] = np.shape(pedestal)
        pedH = pedestal[0, :]
        pedW = pedestal[1, :]
        pedP = pedestal[2, :]
        x = np.linspace(0, col, col)
        ax1.plot(x, pedH)
        ax1.set_title('Pedestal Height')
        ax2.plot(x, pedW)
        ax2.set_title('Pedestal Width')
        ax3.plot(x, pedP)
        ax3.set_title('Pedestal Gradient')
        # plt.pause(0.1)

    def JointData(self,results):
        DataLen = self.MainExc.Counts
        pedestalinfo = np.zeros((3, DataLen))
        for dictK in results.keys():
            clientInfo = results[dictK]
            offset = clientInfo['offset']
            peds = clientInfo['peds']
            # print('peds:',peds)
            if len(peds) > 1:
                pedestalinfo[:, offset[0]:offset[1]] = np.array(peds)
            else:
                pedestalinfo[:, offset[0]:offset[1]] = 0
        return pedestalinfo

    def CauMain(self):
        print('发送任务请求...')
        fileInfo = {}
        fileInfo['FileType'] = 'mat'
        fileInfo['FileData'] = 'dataful86787.mat'
        datas=self.MainExc.SubnitMission(fileInfo)

        JoinClientInfo = json.loads(datas)
        ConsentJoinClient = JoinClientInfo['ConsentJoinClient']
        MissionPerCom = JoinClientInfo['MissionPerCom']
        self.counts = self.MainExc.Counts
        j = 0
        print('总任务数：'+str(MissionPerCom[-1][-1]))
        for i in ConsentJoinClient:
            InfoList = i.split('#')
            IPS = InfoList[0]
            CPUS = int(InfoList[1])
            MissionPer = ((MissionPerCom[j][1]-MissionPerCom[j][0])/MissionPerCom[-1][-1])*100
            print('IP：'+IPS+'\tcpu核数: '+str(CPUS)+'\t任务占比：'+str(MissionPer)+'%')
            j=j+1
        self.CurResult()
        plt.show()

    def FindNodes(self):
        print('寻找计算节点...')
        NodesInfo = self.MainExc.ConsultNodesNums()
        cpuCounts = 0
        for i in NodesInfo:
            info = i.split('#')
            ip = info[0]
            cpus = info[1]
            cpuCounts = cpuCounts+int(cpus)
            print('IP: {} CPU核数: {}'.format(ip,cpus))
        print('共 {} 个节点，CPU核心总数为: {}'.format(str(len(NodesInfo)),str(cpuCounts)))

    def CurResult(self):
        t1 = time.time()
        usingtimes = {}
        showdataJudge = []
        while True:
            curs = self.MainExc.CurResult()
            # print('teeee',curs)
            IP = []
            Progresses = []
            os.system('cls')
            for dictK in curs.keys():
                if dictK not in usingtimes.keys():
                    usingtimes[dictK] = '0'
                clientInfo = curs[dictK]
                ip = clientInfo['ip']
                progress = round(clientInfo['finshRate'], 2)
                if progress != 100:
                    usingtimes[dictK] = str(time.time() - t1)
                else:
                    if ip not in showdataJudge:
                        # datasTemp = JointData(curs)
                        # showResult(datasTemp)
                        showdataJudge.append(ip)
                IP.append(ip)
                Progresses.append(progress)
                offset = clientInfo['offset']
                curcount = str(clientInfo['CurCount'])
                lenoffset = str(offset[1] - offset[0])

                predictUseTime = '000'
                if int(curcount) != 0:
                    pretT = float(usingtimes[dictK])*(int(lenoffset)/int(curcount))
                    predictUseTime = str(round(pretT,2))

                print(
                    'IP: ' + ip + ' 进度：' + str(progress) + ' % ({}/{}) Time: '.format(curcount, lenoffset) + usingtimes[
                        dictK]+' 预计总用时：{}'.format(predictUseTime))
            if sum(Progresses) == len(IP) * 100 and len(Progresses)>0:
                print('Using Time: ' + str(time.time() - t1))
                Pedstals = self.JointData(curs)
                self.SaveDatatoMat(Pedstals)
                self.showResult(Pedstals)
                time.sleep(2)
                break
            time.sleep(20)





Opercont = ''
with open('Oper.txt','r') as f:
    Opercont=f.read()
print(Opercont)
oper = CCOper()
oper.JudgeServerRun()

while True:
    Op = input('->')
    if Op == '1':
        oper.FindNodes()
    elif Op == '2':
        oper.CauMain()


