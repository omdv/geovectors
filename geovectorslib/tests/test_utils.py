# pylint: disable=redefined-outer-name,unused-variable,expression-not-assigned,singleton-comparison

import numpy as np

from geovectorslib import utils


def describe_wrap90deg():
    def convert():
        np.testing.assert_array_equal(
            utils.wrap90deg(np.array([-91, -90, -89, 0, 89, 90, 92])),
            np.array([-89, -90, -89, 0, 89, 90, 88]),
        )


def describe_wrap180deg():
    def convert():
        np.testing.assert_array_equal(
            utils.wrap180deg(np.array([-181, -180, -170, 0, 170, 180, 182])),
            np.array([179, -180, -170, 0, 170, -180, -178]),
        )


def describe_wrap360deg():
    def convert():
        np.testing.assert_array_equal(
            utils.wrap360deg(np.array([-365, -360, -350, -1, 0, 300, 360, 365])),
            np.array([355, 0, 10, 359, 0, 300, 360, 5]),
        )
