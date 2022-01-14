
import socket

class LocalInfo():
    def __init__(self):
        self.ConfigurationFile = 'configuration.txt'
        self.HostIp = self.get_host_ip(self.ConfigurationFile)
        self.HostCpuNums=self.get_CpuNums()
        self.HostName = self.get_host_name()

    def get_configurationInfo(self):
        with open(self.ConfigurationFile,'r') as f:
            Configuration = f.read()
        Configuration = Configuration.split('\n')
        IP = Configuration[0].split('#')[-1]
        JoinCauNet = Configuration[1].split('#')[-1]
        OpenCPUS = Configuration[2].split('#')[-1]
        ServerIP = Configuration[3].split('#')[-1]
        ServerPoint = Configuration[4].split('#')[-1]
        return [IP,JoinCauNet,OpenCPUS,ServerIP,ServerPoint]
    def get_host_ip(sef,ConfigurationFile):

        with open(ConfigurationFile,'r') as f:
            Configuration = f.read()
        Configuration = Configuration.split('\n')
        IP = Configuration[0].split('#')[-1]
        if IP =='0':
            hostname = socket.gethostname()
            # 获取本机IP
            ip = socket.gethostbyname(hostname)
        else:
            ip = IP
        return ip
    def get_CpuNums(self):
        import multiprocessing
        processnums = multiprocessing.cpu_count()
        return processnums
    def get_host_name(self):
        import getpass
        CurHostName=getpass.getuser()
        return CurHostName


# if __name__ == '__main__':
#     print(get_host_ip())
#     import multiprocessing
#     processnums=multiprocessing.cpu_count()
#     print(processnums)