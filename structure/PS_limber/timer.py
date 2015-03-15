'''from timer import Timer
from redis import Redis
rdb = Redis()

with Timer() as t:
    rdb.lpush("foo", "bar")
print "=> elasped lpush: %s s" % t.secs

with Timer() as t:
    rdb.lpop("foo")
print "=> elasped lpop: %s s" % t.secs'''





import time

class Timer(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.verbose:
            print 'elapsed time: %f ms' % self.msecs