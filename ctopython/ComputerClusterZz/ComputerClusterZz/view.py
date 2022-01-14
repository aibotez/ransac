from django.http import HttpResponseRedirect
from django.http import StreamingHttpResponse,FileResponse
from django.http import HttpResponse,JsonResponse
import json,requests
import numpy as np
import threading

import GetLocalInfo
import Main

Server_IP = '222.204.48.43'
Server_Point = '9988'
ClientInfoDocument = 'ClientInfo.txt'
LocalInfos = GetLocalInfo.LocalInfo()
MainExc = Main.ExcuMain()
ClientInfos=MainExc.ClientInfos
# print(ClientInfos)
ClientInfosWaitInsert=[]
if ClientInfos == None:
    ClientInfos = []




# def ReadLocalClientInfo():
#     with open(ClientInfoDocument,'r+',encoding='utf-8') as f:
#         Content = f.read()
#     ClientInfos0 = Content.split('\n')
#     del ClientInfos0[-1]
#     return ClientInfos0

def UpdateClientInfo():
    global ClientInfosWaitInsert,ClientInfos
    while True:
        if len(ClientInfosWaitInsert) >0:
            with open(ClientInfoDocument, 'a+', encoding='utf-8') as f:
                for i in ClientInfosWaitInsert:
                    f.write(i+'\n')
            ClientInfosWaitInsert = []

UpdateClientInfoTread = threading.Thread(target=UpdateClientInfo)
UpdateClientInfoTread.setDaemon(True)
UpdateClientInfoTread.start()



ClientInfoDocumentWrJudge = 0
def GetClientInfos(request):
    requestcontent = request.body
    requestcontent = json.loads(requestcontent)

    hostname = requestcontent['HostName']
    HostNames = [i.split('#')[2] for i in ClientInfos]
    clientinfo = requestcontent['clientinfo']
    # print('ClientInfos:',ClientInfos)
    if clientinfo not in ClientInfos and hostname not in HostNames:
        ClientInfos.append(clientinfo)
        ClientInfosWaitInsert.append(clientinfo)
    return HttpResponse('OK')

def SendJoinMissionRequest(request):
    return HttpResponse('1')

def JudgeServerRun(request):
    return HttpResponse('1')



CurResult = {}
def SubmitMission(request):
    global CurResult
    CurResult = {}
    requestcontent = request.body
    requestcontent = json.loads(requestcontent)
    R0 = np.array(requestcontent['R'])
    ne0 = np.array(requestcontent['ne'])
    [row,col] = np.shape(R0)
    # print('col:',col)
    ConsentJoinClient = MainExc.SendJoinMissionRequest()
    ConsentJoinClientInfo = MainExc.DisMission(ConsentJoinClient,col)
    SplitData = MainExc.SplitData(ConsentJoinClientInfo,R0,ne0)
    MainExc.SendSplitData(SplitData)
    # SendSplitDataTread = threading.Thread(target=MainExc.SendSplitData,args=(SplitData,))
    # SendSplitDataTread.setDaemon(True)
    # SendSplitDataTread.start()
    # MainExc.SendSplitData(SplitData)
    # print(ConsentJoinClientInfo)
    return HttpResponse(json.dumps(ConsentJoinClientInfo), content_type="application/json")


def ConsultNodesNums(request):
    ConsentJoinClient = MainExc.SendJoinMissionRequest()
    return HttpResponse(json.dumps(ConsentJoinClient), content_type="application/json")

def SendSplitData(request):
    requestcontent = request.body
    requestcontent = json.loads(requestcontent)
    MainExc.Compact(requestcontent)
    # CompactTread = threading.Thread(target=MainExc.Compact,args=(requestcontent,))
    # CompactTread.setDaemon(True)
    # CompactTread.start()
    # R=requestcontent['R']
    # ne = requestcontent['ne']
    # MainExc.ExcuteComp(requestcontent)
    return HttpResponse('1')

# def SubmitMissionResult(request):
#     global CurCounts,a
#     requestcontent = request.body
#     requestcontent = json.loads(requestcontent)
#     ai = requestcontent['a']
#     offset = requestcontent['offset']
#     ai = np.array(ai)
#     [row,col] = np.shape(ai)
#     a[:,offset[0]:offset[1]] = ai
#     Counts = MainExc.Counts
#     CurCounts = CurCounts + col
#
#     print(offset,CurCounts)
#     print(a)
#     return HttpResponse('1')



def GetCurResultsfromClient(request):
    global CurResult
    requestcontent = request.body
    # print(requestcontent)
    requestcontent = json.loads(requestcontent)
    for dictK in requestcontent.keys():
        CurResult[dictK] = requestcontent[dictK]
    return HttpResponse('1')


def GetCurResult(request):
    global CurResult
    return HttpResponse(json.dumps(CurResult), content_type="application/json")

def test(request):
    prog = MainExc.prog()
    return HttpResponse(str(prog))

