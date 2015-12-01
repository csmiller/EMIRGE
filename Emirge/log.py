"""Helper functions for log output"""

from functools import wraps
import logging
from time import time
from datetime import timedelta

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(message)s'
)


def timed(msg):
    """
    Decorator adding a log message around a function.
    Once the function has completed, the execution time is shown.
    The first word of msg is repeated in the closing msg.
    """

    first_word = msg.split(" ", 1)[0]
    finished = "DONE " + first_word + "  [%s]"

    def decorator(fn):
        @wraps(fn)
        def inner(*args, **kwargs):
            logging.info(msg + "...")
            start_time = time()
            result = fn(*args, **kwargs)
            end_time = time()
            duration = timedelta(seconds=end_time - start_time)
            logging.info(finished % duration)
            return result
        return inner
    return decorator


error = logging.error
info = logging.info
warning = logging.warning
