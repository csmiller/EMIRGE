"""Helper functions for log output"""

from functools import wraps
import logging
from time import time
from datetime import timedelta
from inspect import getargspec


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(message)s'
)

DEBUG = logging.debug        # just for debug
INFO = logging.info          # just FYI
WARNING = logging.warning    # something to consider
ERROR = logging.error        # something went wrong, can continue
CRITICAL = logging.critical  # something went wrong, need to abort
LOG = logging.log            # giving level

debug = DEBUG
info = INFO
warning = WARNING
error = ERROR
critical = CRITICAL
log = LOG


def setup(quiet=False, debug=False):
    logger = logging.getLogger()
    if debug:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.WARNING)


class Timed(object):
    def __init__(myself, msg, level=logging.INFO, *args, **kwargs):
        myself.msg = msg.format(*args, **kwargs)
        myself.level = level
        myself.start_time = time()
        LOG(level, myself.msg + "...")

    def done(myself):
        LOG(myself.level, "DONE {msg}  [{time}]".format(
            msg=myself.msg.split(" ", 1)[0],
            time=timedelta(seconds=time() - myself.start_time)
        ))

    def failed(myself, msg):
        LOG(logging.ERROR, "ERROR: {msg}  [{time}]".format(
            msg=msg,
            time=timedelta(seconds=time() - myself.start_time)
        ))


def debugTimed(msg, *args, **kwargs):
    return Timed(msg, level=logging.DEBUG, *args, **kwargs)


def infoTimed(msg, *args, **kwargs):
    return Timed(msg, level=logging.INFO, *args, **kwargs)


def warningTimed(msg, *args, **kwargs):
    return Timed(msg, level=logging.WARNING, *args, **kwargs)

DEBUG_timed = debugTimed
INFO_timed = infoTimed
WARNING_timed = warningTimed


def timed(msg, *msg_args):
    """
    Decorator adding a log message around a function.
    Once the function has completed, the execution time is shown.
    The first word of msg is repeated in the closing msg.
    """

    first_word = msg.split(" ", 1)[0]
    finished = "DONE " + first_word + "  [%s]"

    def decorator(fn):
        env = {}
        try:
            argspec = getargspec(fn)
            if argspec.defaults is not None:
                for key, value in zip(reversed(argspec.args),
                                      reversed(argspec.defaults)):
                    env[key] = value
        except TypeError:
            argspec = None

        @wraps(fn)
        def wrapper(*args, **kwargs):
            lenv = dict(env)
            if argspec is not None:
                argnames = argspec.args
            else:
                argnames = ["self", "arg1", "arg2", "arg3"]

            for key, value in zip(argnames, args):
                lenv[key] = value

            try:
                info(msg.format(**lenv))
            except AttributeError as e:
                warning('Failed to format log msg <<<{}>>> with args {}: "{}"'
                     .format(msg, lenv.keys(), e.message))
            start_time = time()
            result = fn(*args, **kwargs)
            end_time = time()
            duration = timedelta(seconds=end_time - start_time)
            logging.info(finished % duration)
            return result
        return wrapper
    return decorator
