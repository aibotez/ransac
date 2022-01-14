import json,requests
import numpy as np
import threading
import GetLocalInfo
import pandas as pd
from scipy import io
import ctopytest

class ExcuMain():
    def __init__(self):
        # self.Server_IP = '222.204.48.43'
        # self.Server_Point = '9988'
        self.LocalInfos = GetLocalInfo.LocalInfo()
        self.ClientInfoDocument = 'ClientInfo.txt'
        self.ClientInfos=self.ReadLocalClientInfo()
        # self.ransac = ctopytest.Test()
        self.Counts = 0
        self.a = None
        self.configurationInfo = self.LocalInfos.get_configurationInfo()
        self.Server_IP = self.configurationInfo[3]
        self.Server_Point = self.configurationInfo[4]

        self.UseCPUs = int(int(self.LocalInfos.HostCpuNums)*int(self.configurationInfo[2])/100)

    def CurResult(self):
        ServerURL = 'http://' + self.Server_IP + ':' + self.Server_Point + '/GetCurResult/'
        res = requests.get(url=ServerURL).text
        Curs = json.loads(res)
        # Curs['progress'] = Curs['CurCounts']/self.Counts
        return Curs

    def ReadLocalClientInfo(self):
        with open(self.ClientInfoDocument, 'r+', encoding='utf-8') as f:
            Content = f.read()
        ClientInfos0 = Content.split('\n')
        del ClientInfos0[-1]
        return ClientInfos0

    # def prog(self):
        # return self.ransac.GetCurProcess()
    def SendCurResult(self,CurPro,peds,fishsatus,offset,finshRate):
        url = 'http://'+self.Server_IP+':'+self.Server_Point+'/SendCurResultstoSer/'
        # print('SendRersultUrl:',url)
        data={}
        # print('LocalInfos.HostIp',self.LocalInfos.HostIp)
        if fishsatus == 'unfish':
            data = {
                'CurCount': CurPro,
                'peds': [],
                'fishsatus': fishsatus,
                'offset': offset,
                'finshRate': finshRate,
                'ip':self.LocalInfos.HostIp,
            }
        elif fishsatus == 'fish':
            data = {
                'CurCount':CurPro,
                'peds':peds,
                'fishsatus':fishsatus,
                'offset':offset,
                'finshRate':finshRate,
                'ip': self.LocalInfos.HostIp,
            }
        DatabyIp = {
            self.LocalInfos.HostIp:data,
        }
        # print('DatabyIp: ',DatabyIp)
        DatabyIp = json.dumps(DatabyIp)
        # print('Data:',DatabyIp)
        res = requests.post(url=url,data=DatabyIp,timeout=3)
    def Compact(self,SplitData):
        # while True:
        #     print('Tr')
        R = SplitData['R']
        ne = SplitData['ne']
        offset = SplitData['offset']
        Rn = np.array(R)
        nen = np.array(ne)
        [row,col] = np.shape(Rn)
        pedsInfo = np.zeros((3, col))

        j = 0
        showpernum = 10
        for i in range(col):
            # print('ii', i)
            Ri = Rn[:,i].tolist()
            nei = nen[:,i].tolist()
            ransac = ctopytest.Test()
            peds=ransac.ExcuteRan(Ri, nei,self.UseCPUs)  # 执行计算
            pedsInfo[:,i] = np.array(peds)
            finshRate = ((i+1)/col)*100
            # self.SendCurResult(i + 1, peds, 'unfish', offset, finshRate)
            if i % showpernum ==0:
                try:
                    self.SendCurResult(i + 1, peds, 'unfish', offset, finshRate)
                except:
                    print('SendCurResult ERROE')
        pedsInfoL = pedsInfo.tolist()
        finshRate = 100
        while True:
            try:
                self.SendCurResult(col, pedsInfoL, 'fish', offset,finshRate)
                break
            except:
                pass
        #
        #
        # self.a = self.ransac.ExcuteRan(R, ne)#执行计算
        # ServerURL = 'http://' + self.Server_IP + ':' + self.Server_Point + '/SubmitMissionResult/'
        # data={
        #     'a':self.a,
        #     'offset':SplitData['offset'],
        # }
        # data = json.dumps(data)
        # res = requests.post(url=ServerURL, data=data).text

    def ExcuteComp(self,SplitData):
        CompTread = threading.Thread(target=self.Compact,args=(SplitData,))
        CompTread.setDaemon(True)
        CompTread.start()
    def SendSplitData(self,SplitData):
        for i in SplitData:
            ip = i['ip']
            R = i['R'].tolist()
            ne = i['ne'].tolist()
            data = {
                'R':R,
                'ne':ne,
                'offset':i['offset'],
            }
            data = json.dumps(data)
            URL = 'http://'+ip + '/SendSplitData/'
            # print(URL)
            try:
                res = requests.get(url=URL,data=data,timeout=3).text
            except:
                pass

    def SplitData(self,ConsentJoinClientInfo,R,ne):
        MissionPerCom = ConsentJoinClientInfo['MissionPerCom']
        ClientInfo = ConsentJoinClientInfo['ConsentJoinClient']

        SpDataInfo = []
        j=0
        for i in MissionPerCom:
            MissionInfoPerCom = {'ip': '', 'R': None, 'ne': None, 'offset': []}
            R1 = R[:,i[0]:i[1]]
            ne1 = ne[:, i[0]:i[1]]
            ip = ClientInfo[j].split('#')[0]
            # print('ippp:',ip)
            MissionInfoPerCom['ip'] = ip
            MissionInfoPerCom['R'] = R1
            MissionInfoPerCom['ne'] = ne1
            MissionInfoPerCom['offset'] = i
            SpDataInfo.append(MissionInfoPerCom)
            j=j+1
        return SpDataInfo

    def DisMission(self,ConsentJoinClient,col):
        self.Counts = col
        IPS = []
        CPUS = []
        for i in ConsentJoinClient:
            InfoList = i.split('#')
            IPS.append(InfoList[0])
            CPUS.append(int(InfoList[1]))
        CPUSacount = sum(CPUS)
        MissionPerCom = []
        CountPer = col/CPUSacount

        OffsetCot = 0
        for i in range(len(CPUS)):
            Cots = int(CountPer*CPUS[i])
            MissionPerCom.append([OffsetCot,OffsetCot+Cots])
            OffsetCot = OffsetCot+Cots+0
        MissionPerCom[-1][-1] = col
        return {'ConsentJoinClient':ConsentJoinClient,'MissionPerCom':MissionPerCom}

    def ConsultNodesNums(self):
        url = 'http://' + self.Server_IP + ':' + self.Server_Point + '/ConsultNodesNums/'
        res = requests.get(url=url).text
        data = json.loads(res)
        return data

    def SendJoinMissionRequest(self):
        ConsentJoinClient = []
        self.ReadLocalClientInfo()
        for i in self.ClientInfos:
            ip = i.split('#')[0]
            URL = 'http://'+ip + '/SendJoinMissionRequest/'
            try:
                res = requests.get(url=URL,timeout=5).text
                if res == '1':
                    ConsentJoinClient.append(i)
            except:
                pass
        return ConsentJoinClient

    def SubnitMission(self,fileInfo):
        ServerURL = 'http://' + self.Server_IP + ':' + self.Server_Point + '/SubmitMission/'

        filetype = fileInfo['FileType']
        if filetype == 'txt':
            fileR = fileInfo['FileR']
            filene = fileInfo['Filene']
            Rn = pd.read_table(fileR, header=None)
            R = np.array(Rn)
            [row, self.Counts] = np.shape(R)
            R = R.tolist()
            nen = pd.read_table(filene, header=None)
            ne = np.array(nen).tolist()
        elif filetype == 'mat':
            filedata = fileInfo['FileData']
            data = io.loadmat(filedata)
            R = data['R']
            try:
                ne = data['ne']
            except:
                ne = data['N']
            [row, self.Counts] = np.shape(R)
            R = R.tolist()
            ne = ne.tolist()
        data={
            'R':R,
            'ne':ne,
        }

        data = json.dumps(data)
        # print(len(data))
        res = requests.post(url=ServerURL, data=data).text
        return res

    def SubmitComputerInfos(self):
        ip = self.LocalInfos.HostIp
        CompCpus = self.LocalInfos.HostCpuNums
        HostName = self.LocalInfos.HostName
        if ':' in ip:
            hostip = ip
        else:
            hostip = ip + ':{}'.format(str(self.Server_Point))
        ComputerInfos = {
            'HostIp': hostip,
            'CPUs': CompCpus,
            'HostName': HostName,
            'clientinfo':hostip+'#'+str(CompCpus)+'#'+HostName,
        }
        ServerURL = 'http://'+self.Server_IP + ':' + self.Server_Point + '/GetClientInfos/'
        data = json.dumps(ComputerInfos)
        res = requests.post(url=ServerURL, data=data).text

        return res