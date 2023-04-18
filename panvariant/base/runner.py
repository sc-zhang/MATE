from os import popen
from panvariant.io.message import Message as Msg


class Runner:
    def __init__(self):
        self.__cmd = ""
        self.__res = ""

    def set_cmd(self, cmd):
        self.__cmd = cmd

    def run(self):
        try:
            with popen(self.__cmd, 'r') as fin:
                res = []
                for _ in fin:
                    res.append(_)
            self.__res = '\n'.join(res)
            return 1
        except Exception as e:
            Msg.error("\tRunning command {} failed: {}".format(self.__cmd, repr(e)))

    def get_result(self):
        return self.__res
