import zmq
import logging
import time
import json
import numpy as np
import time
from copy import deepcopy
from functools import reduce

VERBOSE=True

try:
    from icecream2 import ic
except ModuleNotFoundError as e:
    pass

def _ic(m):
    if not VERBOSE:
        return
    try:
        ic(m)
    except NameError as e:
        print(m)

def cache(caching, threshold = 1.0e-6):

    # caching should be list of strings:

    if not type(caching) in (list, tuple):
        raise TypeError("Variable caching has to be list or tuple")

    for t in caching:
        if not type(t) == str:
            raise TypeError("Variable caching has to be list or tuple of str")

    # Why int that is crazy!
    if not type(threshold) in (int, float):
        raise TypeError("Variable threshold has to be int or float")

    def wrapper(f):
        def func(self, data):
            if not self.last:
                self.eval(data)
            else:
                if reduce(lambda x, y: x or y,
                          [(np.array(data[val]) -  np.array(self.last[val])).any() > threshold for val in caching]):
                    self.eval(data)
            return f(self, data)
        return func
    return wrapper

class Jastrow():

    def __init__(self):
        self.last = None

    def eval(self, data):
        self._value = list(np.zeros(len(data["elec_coord"])//3))
        self._grad = list(np.zeros(len(data["elec_coord"])))
        self._laplace = list(np.zeros(len(data["elec_coord"])//3))
        self.last = deepcopy(data)

    @cache(caching = ["elec_coord"])
    def value(self, data):
        return self._value

    @cache(caching = ["elec_coord"])
    def grad(self, data):
        return self._grad

    @cache(caching = ["elec_coord"])
    def laplace(self, data):
        return self._laplace


class Server():

    def __init__(self, jastrow, port=5555):
        self.s = []
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.REP)
        self.socket.setsockopt(zmq.RCVBUF, 1001)
        self.socket.setsockopt(zmq.RCVHWM, 1001)
        self.socket.bind(f"tcp://*:{port}")
        self.running = True
        self.jastrow = jastrow
        self.count = 0

    def run(self):
        while self.running:
            message = self.socket.recv().decode()
            #time.sleep(0.3)
            if message.strip().lower() == "stop":
                self.running = False
                ret = f"Exiting"
                self.socket.send(ret.encode())
                self.socket.close()
                break
            try:
                _ic(message)
                data = json.loads(message)
                _ic(f"request {self.count}")
                if data["type"] == "value":
                    message = " ".join([str(x) for x in self.jastrow.value(data) ])
                elif data["type"] == "grad":
                    message = " ".join([str(x) for x in self.jastrow.grad(data) ])
                elif data["type"] == "laplace":
                    message = " ".join([str(x) for x in self.jastrow.laplace(data) ])
                else:
                    message("Unknown type")
            except Exception as e:
                _ic(e)
            ret = message + " "*(100-len(message))
            _ic(message.strip())
            self.socket.send(ret.encode())
            _ic("message sended")
            self.count += 1
        self.socket.close()
        self.finalize()

    def finalize(self):
        for s in self.s:
            s.join()

if __name__ == "__main__":
    s = Server(Jastrow(), port = 6543)
    s.run()

