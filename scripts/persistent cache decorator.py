# -*- coding: utf-8 -*-
from functools import wraps


def memo(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        x = f(*args, **kwds)
        # make decisions here on when to retrieve cache depending on env variables etc..
        # right now we only have a saved value for 1+2+3
        #
        if x.key == "1+2+3":  # we know 1+2+3 is 6, we can use the cached value
            print("cached value returned")
            return 6
        return x.sum()  # <-- this method has to have the same name in all

    return wrapper


class Example:
    """Docstring for Example class"""

    @memo
    def __init__(self, a=1, b=2, c=3):
        self.a = a
        self.b = b
        self.c = c

    def sum(self):
        print("value calculated")
        return self.a + self.b + self.c


print(Example(1, 2, 3))

print(Example(1, 2, 4))

print(Example.__doc__)
