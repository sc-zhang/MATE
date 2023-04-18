from os import popen


class DepCheck:
    def __init__(self):
        pass

    @staticmethod
    def check(cmd):
        res = []
        with popen(cmd, 'r') as fin:
            for line in fin:
                res.append(line.strip())

        if res:
            return True
        return False
