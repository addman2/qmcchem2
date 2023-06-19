from Jastrow_evaluator import Jastrow

class JastrowQMCkl(Jastrow):

    def __init__(self, trexio_file, *args, **kwargs):
        import trexio
        super().__init__(self, *args, **kwargs)

        self.trexio_file = trexio_file
        self.qmckl_instance = None

    def eval(self, data):
        # self.value = self.qmckl_instance.get_val(...)
        self.value = list(np.zeros(len(data["nucl_coord"])//3))
        self.grad = list(np.zeros(len(data["nucl_coord"])))
        self.laplace = list(np.zeros(len(data["nucl_coord"])//3))
        self.last = data

