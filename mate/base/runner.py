from subprocess import Popen, PIPE
from mate.io.message import Message as Msg


class Runner:
    def __init__(self):
        self.__cmd = ""
        self.__res = ""
        self.__err = ""

    def set_cmd(self, cmd):
        self.__cmd = cmd

    def run(self):
        try:
            p = Popen(self.__cmd, stdout=PIPE, stderr=PIPE, shell=True, encoding='utf-8')
            self.__res, self.__err = p.communicate()
            return 1
        except Exception as e:
            Msg.error("\tRunning command {} failed: {}".format(self.__cmd, repr(e)))

    def get_result(self):
        return self.__res

    def get_err(self):
        return self.__err
