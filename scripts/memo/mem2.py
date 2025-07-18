# -*- coding: utf-8 -*-
from functools import partial


class memoize(object):
    def __init__(self, fn):
        self.fn = fn

    def __call__(self, func):
        def wrappee(*args, **kwargs):
            obj = args[0]
            key = obj.n  # key depends on obj proerties
            print(args)
            cache = {3: 42, 2: 55}
            try:
                res = cache[key]
            except KeyError:
                res = cache[key] = func(*args, **kwargs)
            return res

        return wrappee


if __name__ == "__main__":

    class doubleclass(object):
        def __init__(self, n):
            self.n = n

        @memoize("hej")
        def res(self):
            return self.n * 2

    x = doubleclass(2)
    print(x.res())
