import zmq
import json
from icecream import ic
import numpy as np

def cache(f):

    def func(self, data):
        if not self.last:
            self.eval(data)
        else:
            if (np.array(data["nucl_coord"]) -  np.array(self.last["nucl_coord"])).any() > 1.0e-6 or \
               (np.array(data["elec_coord"]) -  np.array(self.last["elec_coord"])).any() > 1.0e-6:
                self.eval(data)
        return f(self, data)
    return func

class Jastrow():

    def __init__(self):
        self.last = None

    def eval(self, data):
        self.value = list(np.zeros(len(data["nucl_coord"])//3))
        self.grad = list(np.zeros(len(data["nucl_coord"])))
        self.laplace = list(np.zeros(len(data["nucl_coord"])//3))
        self.last = data

    @cache
    def value(self, data):
        return self.value

    @cache
    def grad(self, data):
        return self.grad

    @cache
    def laplace(self, data):
        return self.laplace

class Server():

    def __init__(self, jastrow):
        self.s = []
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.REP)
        self.socket.bind("tcp://*:5555")
        self.running = True
        self.jastrow = jastrow

    def run(self):
        while self.running:
            message = self.socket.recv().decode()
            if message.strip().lower() == "stop":
                self.running = False
                ret = f"Exiting"
                self.socket.send(ret.encode())
                self.socket.close()
                break
            try:
                ic(message)
                data = json.loads(message)
                ic(data)
                if data["type"] == "value":
                    message = " ".join([str(x) for x in self.jastrow.value(data) ])
                if data["type"] == "grad":
                    message = " ".join([str(x) for x in self.jastrow.grad(data) ])
                if data["type"] == "laplace":
                    message = " ".join([str(x) for x in self.jastrow.laplace(data) ])
            except Exception as e:
                ic(e)
            ret = message
            self.socket.send(ret.encode())
        self.socket.close()
        self.finalize()

    def finalize(self):
        for s in self.s:
            s.join()

if __name__ == "__main__":
    s = Server(Jastrow())
    s.run()
