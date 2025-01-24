# -*- coding: utf-8 -*-
import time, pickle, tempfile, os


class AutoPickleMeta(type):
    def __call__(cls, *args):
        t = (int(time.time()) // 10000) * 10000
        h = hash(args)
        fn = "%s/%s-%i-%i.pickle" % (tempfile.gettempdir(), cls.__name__, t, h)

        if os.path.exists(fn):
            # File exists, so load the cPickle and return
            with open(fn, "rb") as f:
                try:
                    return pickle.load(f)

                except pickle.UnpicklingError:
                    # If error occurs assume cPickle file is corrupt,
                    # and create a new object
                    f.close()
                    return cls._do_pickle(fn, args)

                except EOFError:
                    # File appears empty, return a new object
                    f.close()
                    return cls._do_pickle(fn, args)
        else:
            return cls._do_pickle(fn, args)

    def _do_pickle(cls, fn, args):
        # Create object, and return
        o = object.__new__(cls, *args)
        o.__init__(*args)
        with open(fn, "wb") as f:
            try:
                pickle.dump(o, f, pickle.HIGHEST_PROTOCOL)

            except:
                # If anything went wrong, delete the pickle file and re-raise
                # exception
                f.close()
                os.remove(fn)
                raise
        return o


class AutoPickle(object):
    __metaclass__ = AutoPickleMeta
