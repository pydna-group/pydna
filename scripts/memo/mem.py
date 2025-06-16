# -*- coding: utf-8 -*-
from functools import partial


class memoize(object):
    def __init__(self, func):
        self.func = func

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self.func
        return partial(self, obj)

    def __call__(self, *args, **kw):
        obj = args[0]
        key = obj.n
        cache = {2: 42}
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kw)
        return res


if __name__ == "__main__":

    class doubleclass(object):
        def __init__(self, n):
            self.n = n

        @memoize
        def res(self):
            return self.n * 2

    x = doubleclass(2)
    print(x.res())
