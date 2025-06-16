# -*- coding: utf-8 -*-
import inspect


def get_id_tuple(
    f, args, kwargs, mark=object()
):  # Some quick'n'dirty way to generate a unique key for an specific call.
    l = [id(f)]
    for arg in args:
        l.append(id(arg))
    l.append(id(mark))
    for k, v in kwargs:
        l.append(k)
        l.append(id(v))
    return tuple(l)


_memoized = {}


def memoize(ex):  # Some basic memoizer
    def decorator(f):
        # print(inspect.isclass(f), inspect.ismethod(f), inspect.isfunction(f))
        def memoized(*args, **kwargs):
            key = get_id_tuple(f, args, kwargs)
            if key not in _memoized:
                _memoized[key] = f(*args, **kwargs)
            import pickle

            x = f(*args, **kwargs)
            print(x.__name__)
            pickle.dumps()
            return _memoized[key]

        return memoized

    return decorator


@memoize("bajs")
class Test(object):
    def __init__(self, somevalue):
        self.somevalue = somevalue


@memoize("phlegm")
def myfunc(x):
    return 2 * x


class Test2(object):
    def __init__(self, somevalue):
        self.somevalue = somevalue

    @memoize("piss")
    def do(self):
        return 2 * self.somevalue


tests = [Test(1), Test(2), Test(2)]
for test in tests:
    print(test.somevalue, id(test))
print()
tests = [myfunc(1), myfunc(2), myfunc(2)]
for test in tests:
    print(test, id(test))
print()
tests = [Test2(1), Test2(2), Test2(2)]
for test in tests:
    print(test.do(), id(test))
