import os,time,requests,json
import threading
import Main
import numpy as np
import pandas as pd

MainExc = Main.ExcuMain()
def ReadDataFromtxt(fileR,filene):
    R0 = pd.read_table(fileR, header=None)
    ne0 = pd.read_table(filene, header=None)
    return [R0,ne0]

def StartProject():
    print('Start..')
    os.system('startproject.cmd')
    # os.system('python manage.py runserver 0.0.0.0:9988')

def JudgeServerRun():
    configurationInfo = MainExc.configurationInfo
    if configurationInfo[1] =='1':
        res = MainExc.SubmitComputerInfos()

        try:
            res = requests.get('http://127.0.0.1:9988/JudgeServerRun').text
        except:
            os.system('start startproject.cmd')
            time.sleep(5)
    else:
        pass
# StartProjectTread = threading.Thread(target=StartProject)
# StartProjectTread.setDaemon(True)
# StartProjectTread.start()


JudgeServerRun()



def CurResult():
    t1 = time.time()
    while True:
        curs=MainExc.CurResult()
        IP = []
        Progresses = []
        os.system('cls')
        usingtimes = {}
        for dictK in curs.keys():
            if dictK not in usingtimes.keys():
                usingtimes[dictK] = '0'
            clientInfo = curs[dictK]
            ip = clientInfo['ip']
            progress = clientInfo['finshRate']
            if progress != 100:
                usingtimes[dictK] = str(time.time() - t1)
            IP.append(ip)
            Progresses.append(progress)
            offset = clientInfo['offset']
            curcount = str(clientInfo['CurCount'])
            lenoffset = str(offset[1]-offset[0])
            print('IP: ' + ip + ' 进度：' + str(progress) + ' ({}/{}) Time: '.format(curcount,lenoffset) + usingtimes[dictK]+'\n')
        if sum(Progresses) == len(IP)*100:
            print('Using Time: ' + str(time.time() - t1))
            break


        # progress = curs['progress']*100
        # print('进度：'+str(progress)+'\tTime: '+str(time.time()-t1))
        # if progress >=100:
        #     break
        time.sleep(5)
Opercont = ''
with open('Oper.txt','r') as f:
    Opercont=f.read()

# os.system('start startproject.cmd')
# time.sleep(2)

print(Opercont)

fileR = 'R0.txt'
filene = 'ne0.txt'
while True:
    Op = input('->')
    if Op == '1':
        print('寻找计算节点...')
        NodesInfo = MainExc.ConsultNodesNums()
        cpuCounts = 0
        for i in NodesInfo:
            info = i.split('#')
            ip = info[0]
            cpus = info[1]
            cpuCounts = cpuCounts+int(cpus)
            print('IP: {} CPU核数: {}'.format(ip,cpus))
        print('共 {} 个节点，CPU核心总数为: {}'.format(str(len(NodesInfo)),str(cpuCounts)))
    elif Op == '2':
        print('发送任务请求...')
        datas=MainExc.SubnitMission(fileR,filene)
        JoinClientInfo = json.loads(datas)
        ConsentJoinClient = JoinClientInfo['ConsentJoinClient']
        MissionPerCom = JoinClientInfo['MissionPerCom']
        j = 0
        print('总任务数：'+str(MissionPerCom[-1][-1]))
        for i in ConsentJoinClient:
            InfoList = i.split('#')
            IPS = InfoList[0]
            CPUS = int(InfoList[1])
            MissionPer = ((MissionPerCom[j][1]-MissionPerCom[j][0])/MissionPerCom[-1][-1])*100
            print('IP：'+IPS+'\tcpu核数: '+str(CPUS)+'\t任务占比：'+str(MissionPer)+'%')
            j=j+1
        CurResult()


