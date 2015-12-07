"""Test Emirge.log"""

from Emirge import log
import logging
from nose.tools import assert_equal


class MockLogHandler(logging.Handler):
    """Mock handler for logging
    Captures log entries so that they can be inspected by unit tests
    """
    def __init__(self, *args, **kwargs):
        self.messages = {
            'debug': [], 'info': [], 'warning': [],
            'error': [], 'critical': []
        }
        super(MockLogHandler, self).__init__(*args, **kwargs)

    def emit(self, rec):
        self.acquire()
        try:
            self.messages[rec.levelname.lower()].append(rec.getMessage())
        finally:
            self.release()

    def reset(self):
        self.acquire()
        try:
            self.messages = {}.fromkeys(self.messages, [])
        finally:
            self.release()


class UseMockLog(object):
    def __init__(self):
        self.logs = MockLogHandler(level='DEBUG')
        logger = logging.getLogger()
        logger.addHandler(self.logs)

    @classmethod
    def setup_class(klass):
        """run once before all tests"""

    def setUp(self):
        """run before each test"""
        self.logs.reset()

    def assert_msg_equal(self, msgno, msg):
        assert_equal(self.logs.messages['info'][msgno], msg)

    def assert_msg_startswith(self, msgno, msg):
        assert_equal(self.logs.messages['info'][msgno][:len(msg)], msg)


class Test_Decorator_timed(UseMockLog):
    prop1 = "default_prop1"
    prop2 = "default_prop2"

    @log.timed("simple msg")
    def simple_msg(self, arg1, arg2="default_arg2"):
        pass

    def test_simple_msg(self):
        self.simple_msg(arg1="arg1set")
        self.assert_msg_equal(0, "simple msg")
        self.assert_msg_startswith(1, 'DONE simple  [0:00:00')

    @log.timed("prop msg {self.prop1}")
    def prop_msg(self,  arg1, arg2="default_arg2"):
        pass

    def test_prop_msg(self):
        self.prop_msg(arg1="arg1set")
        self.assert_msg_equal(0, 'prop msg ' + self.prop1)
        self.assert_msg_startswith(1, 'DONE prop  [0:00:00')

    @log.timed("prop2 msg {self.prop2} / {self.prop1}")
    def prop2_msg(self):
        pass

    def test_prop2_msg(self):
        self.prop2_msg()
        self.assert_msg_equal(0,
                              'prop2 msg ' + self.prop2 + ' / ' + self.prop1)
        self.assert_msg_startswith(1, 'DONE prop2  [0:00:00')


class Test_Class_Timed(UseMockLog):
    prop1 = "prop1default"

    def test_simple_msg(self):
        timed_section = log.infoTimed("Simple msg")
        timed_section.done()
        self.assert_msg_equal(0, 'Simple msg...')
        self.assert_msg_startswith(1, 'DONE Simple  [0:00:00')

    def test_prop_msg(self):
        timed_section = log.Timed("Prop msg {self.prop1}", self=self)
        timed_section.done()
        self.assert_msg_equal(0, 'Prop msg ' + self.prop1 + "...")
        self.assert_msg_startswith(1, 'DONE Prop  [0:00:00')
