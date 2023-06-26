import sys
import pytest
from copy import deepcopy

sys.path.append('..')
from Jastrow_evaluator import *

@pytest.fixture
def data():

    yield { "elec_coord": [0.0, 0.0, 0.0],
            "nucl_coord": [0.0, 0.0, 0.0],
            "type": "value" }

def test_cache_BadTreshold1():

    with pytest.raises(TypeError):
        class A():

            def __init__(self, *args, **kwargs):
                pass

            @cache(threshold = "a")
            def eval(data):
                pass

def test_cache_BadTreshold2():

    class A():

        def __init__(self, *args, **kwargs):
            pass

        @cache(caching = ["x"], threshold = 1)
        def eval(data):
            pass

def test_cache_BadTreshold3():

    class A():

        def __init__(self, *args, **kwargs):
            pass

        @cache(caching = ["x"], threshold = 1.0)
        def eval(data):
            pass

def test_cache_BadCache1():

    with pytest.raises(TypeError):
        class A():

            def __init__(self, *args, **kwargs):
                pass

            @cache
            def eval(data):
                pass

def test_cache_BadCache2():

    with pytest.raises(TypeError):
        class A():

            def __init__(self, *args, **kwargs):
                pass

            @cache(caching = 1)
            def eval(data):
                pass

def test_cache_BadCache3():

    with pytest.raises(TypeError):
        class A():

            def __init__(self, *args, **kwargs):
                pass

            @cache(caching = [1])
            def eval(data):
                pass

def test_cache_BadCache4():

    with pytest.raises(TypeError):
        class A():

            def __init__(self, *args, **kwargs):
                pass

            @cache(caching = ["nucl_coord", 1])
            def eval(data):
                pass

def test_cache_GoodCache1():

    class A():

        def __init__(self, *args, **kwargs):
            pass

        @cache(caching = ["nucl_coord"])
        def eval(data):
            pass

def test_cache_GoodCache2():

    class A():

        def __init__(self, *args, **kwargs):
            pass

        @cache(caching = ["elec_coord", "nucl_coord"])
        def eval(data):
            pass

def test_cache_Return1(data):

    class A():

        def __init__(self, *args, **kwargs):
            self.last = None
            self.count = 0

        def eval(self, data):
            self.count += 1
            self.last = deepcopy(data)
            raise Exception()

        @cache(caching = ["elec_coord", "nucl_coord"])
        def value(self, data):
            return self.value

    a = A()

    with pytest.raises(Exception):
        a.value(data)

    assert a.count == 1

    a.value(data)

    assert a.count == 1

def test_cache_Return2(data):

    class A():

        def __init__(self, *args, **kwargs):
            self.last = None
            self.count = 0

        def eval(self, data):
            self.count += 1
            self.last = deepcopy(data)
            raise Exception()

        @cache(caching = ["elec_coord", "nucl_coord"])
        def value(self, data):
            return self.value

    a = A()
    assert a.last is None
    print(a.last)

    with pytest.raises(Exception):
        a.value(data)

    assert a.count == 1

    data["elec_coord"] = [0.0, 0.0, 1.0]

    with pytest.raises(Exception):
        a.value(data)

    assert a.count == 2

def test_cache_Return3(data):

    class A():

        def __init__(self, *args, **kwargs):
            self.last = None
            self.count = 0

        def eval(self, data):
            self.count += 1
            self.last = deepcopy(data)
            raise Exception()

        @cache(caching = ["nucl_coord"])
        def value(self, data):
            return self.value

    a = A()
    assert a.last is None
    print(a.last)

    with pytest.raises(Exception):
        a.value(data)

    assert a.count == 1

    data["elec_coord"] = [0.0, 0.0, 1.0]

    a.value(data)

    assert a.count == 1
