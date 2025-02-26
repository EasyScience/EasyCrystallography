# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

import pytest

from easyscience import global_object
from easyscience.models.polynomial import Line
from easyscience.Objects.Groups import BaseCollection
from easyscience.Objects.variable import DescriptorStr
from easyscience.Objects.variable import Parameter
from easycrystallography.Components.Site import Atoms, Site, _SITE_DETAILS
from easycrystallography.io.star_base import ItemHolder
from easycrystallography.io.star_base import StarLoop
from easycrystallography.io.star_base import StarSection
import gc


@pytest.mark.parametrize(
    "value, variance, precision, expected",
    (
        (1.234560e05, 1.230000e02, 1, "123460(10)"),
        (1.234567e01, 1.230000e-03, 2, "12.346(35)"),
        (1.234560e-01, 1.230000e-04, 3, "0.1235(111)"),
        (1.234560e-03, 1.234500e-08, 4, "0.0012346(1111)"),
        (1.234560e-05, 1.234000e-07, 1, "0.0000(4)"),
    ),
    ids=[
        "1.234560e+05 +- 1.230000e+02 @1",
        "1.234567e+01 +- 1.230000e-03 @2",
        "1.234560e-01 +- 1.230000e-04 @3",
        "1.234560e-03 +- 1.234500e-08 @4",
        "1.234560e-05 +- 1.234000e-07 @1",
    ],
)
def test_ItemHolder_with_error(value, variance, precision, expected):
    p = Parameter("p", value, variance=variance)
    s = ItemHolder(p, decimal_places=precision)
    assert str(s) == expected


@pytest.mark.parametrize("fixed", (True, False), ids=["fixed", "not fixed"])
@pytest.mark.parametrize(
    "value, precision, expected",
    (
        (1.234560e05, 1, "123456.0"),
        (1.234567e01, 2, "12.35"),
        (1.234560e-01, 3, "0.123"),
        (1.234560e-03, 4, "0.0012"),
        (1.234560e-05, 1, "0.0"),
        (1.234560e-05, 5, "0.00001"),
    ),
    ids=[
        "1.234560e+05 @1",
        "1.234567e+01 @2",
        "1.234560e-01 @3",
        "1.234560e-03 @4",
        "1.234560e-05 @1",
        "1.234560e-05 @5",
    ],
)
# FAILED tests/unit_tests/io/test_star.py::test_ItemHolder_fixed[1.234560e-01 @3-fixed] - RuntimeError: dictionary changed size during iteration
# Adding map clear because of these errors happening ONLY in 3.12
def test_ItemHolder_fixed(fixed, value, precision, expected):
    global_object.map._clear()
    p = Parameter("p", value, fixed=fixed)
    s = ItemHolder(p, decimal_places=precision)
    if not p.fixed:
        expected += "()"
    assert str(s) == expected


@pytest.mark.parametrize("cls", [DescriptorStr])
def test_ItemHolder_str(cls):
    v = cls("v", "fooooooooo")
    s = ItemHolder(v)
    assert str(s) == "fooooooooo"


def test_StarSection():
    l = Line(2, 3)
    s = StarSection(l)
    expected = "_m   2.00000000()\n_c   3.00000000()\n"
    assert str(s) == expected


def test_StarLoop():
    gc.collect()
    l1 = Line(2, 3)
    l2 = Line(4, 5)

    ps = BaseCollection("LineCollection", l1, l2)
    s = StarLoop(ps)

    expected = (
        "loop_\n _m\n _c\n  2.00000000()  3.00000000()\n  4.00000000()  5.00000000()"
    )

    assert str(s) == expected
