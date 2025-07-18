# -*- coding: utf-8 -*-
def deco(func):
    def __init__(self, bird):
        self.bird = bird

    def inner(*args):
        print("DECORATED: args={}".format(args))
        func(*args)

    return inner


class Class:
    def __init__(self):
        pass

    @deco("and")
    def method(self, param):
        print("PARAM is {}".format(param))


@deco("anka")
def funk(a, b, c):
    print("{} {} {}".format(a, b, c))


x = Class()

x.method("X")

funk(1, 2, 3)
