
import socket


class LocalInfo():
    def __init__(self):
        self.HostIp = self.get_host_ip()
        self.HostCpuNums=self.get_CpuNums()
    def get_host_ip(sef):
        """
        查询本机ip地址
        :return:
        """
        # 获取计算机名称
        hostname = socket.gethostname()
        # 获取本机IP
        ip = socket.gethostbyname(hostname)
        return ip
    def get_CpuNums(self):
        import multiprocessing
        processnums = multiprocessing.cpu_count()
        return processnums
# if __name__ == '__main__':
#     print(get_host_ip())
#     import multiprocessing
#     processnums=multiprocessing.cpu_count()
#     print(processnums)