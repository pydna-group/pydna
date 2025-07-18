# -*- coding: utf-8 -*-
from functools import wraps


def with_print(func):
    """Decorate a function to print its arguments."""

    @wraps(func)
    def my_func(*args, **kwargs):
        print(args, kwargs)
        import pickle

        pickle.dumps(func(*args, **kwargs))
        return func(*args, **kwargs)

    return my_func


@with_print
class Test(object):
    def __init__(self, somevalue):
        self.somevalue = somevalue


a = Test(123)
