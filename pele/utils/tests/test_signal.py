from __future__ import print_function
import unittest

from pele.utils.events import Signal


class Model(object):
    def __init__(self, value):
        self.__value = value
        self.changed = Signal()

    def set_value(self, value):
        self.__value = value
        self.changed()  # Emit signal

    def get_value(self):
        return self.__value


class View(object):
    def __init__(self, model):
        self.model = model
        model.changed.connect(self.model_changed)
        self.times_changed = 0

    def model_changed(self):
        self.times_changed += 1
        print(("   New value:", self.model.get_value()))


class TestSignal(unittest.TestCase):
    def test1(self):
        print("Beginning Tests:")
        model = Model(10)
        view1 = View(model)
        view2 = View(model)
        view3 = View(model)

        print("Setting value to 20...")
        model.set_value(20)
        self.assertEqual(view1.times_changed, 1)

        print("Deleting a view, and setting value to 30...")
        del view1
        model.set_value(30)
        self.assertEqual(view2.times_changed, 2)

        print("Clearing all listeners, and setting value to 40...")
        model.changed.clear()
        model.set_value(40)
        self.assertEqual(view2.times_changed, 2)

        print("Testing non-member function...")

        self.bar_called = False

        def bar():
            print("   Calling Non Class Function!")
            self.bar_called = True

        model.changed.connect(bar)
        model.set_value(50)
        self.assertEqual(view2.times_changed, 2)
        self.assertEqual(view3.times_changed, 2)
        self.assertTrue(self.bar_called)


if __name__ == "__main__":
    unittest.main()

