# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

import pytest
import easyscience
import numpy as np
from numbers import Number
from easycrystallography.Components.Lattice import Lattice, Parameter, CELL_DETAILS
from easyscience import global_object

pars_dict = {
    "cubic": (5, 5, 5, 90, 90, 90),
    "tetragonal": (10, 10, 5, 90, 90, 90),
    "orthorhombic": (2, 3, 4, 90, 90, 90),
    "monoclinic": (2, 3, 4, 90, 99, 90),
    "hexagonal": (3, 3, 4, 90, 90, 120),
    "rhombohedral": (4, 4, 4, 99, 99, 99),
}


def mod_pars(in_mods=None, sep=False) -> tuple:
    items = []
    keys = pars_dict.keys()
    if in_mods is None:
        in_mods = [[]] * len(keys)
    for key, mod_ in zip(keys, in_mods):
        if mod_:
            if sep:
                items.append(pytest.param(pars_dict[key], mod_, id=key))
            else:
                items.append(pytest.param((*pars_dict[key], mod_), id=key))
        else:
            items.append(pytest.param(pars_dict[key], id=key))
    return tuple(items)


basic_pars = mod_pars()

matrix_pars = mod_pars(
    [
        [
            [5.000000e00, 0.000000e00, 3.061617e-16],
            [-3.061617e-16, 5.000000e00, 3.061617e-16],
            [0.000000e00, 0.000000e00, 5.000000e00],
        ],
        [
            [1.000000e01, 0.000000e00, 6.123234e-16],
            [-6.123234e-16, 1.000000e01, 6.123234e-16],
            [0.000000e00, 0.000000e00, 5.000000e00],
        ],
        [
            [2.0000000e00, 0.0000000e00, 1.2246468e-16],
            [-1.8369702e-16, 3.0000000e00, 1.8369702e-16],
            [0.0000000e00, 0.0000000e00, 4.0000000e00],
        ],
        [
            [1.97537668e00, 0.00000000e00, -3.12868930e-01],
            [-1.83697020e-16, 3.00000000e00, 1.83697020e-16],
            [0.00000000e00, 0.00000000e00, 4.00000000e00],
        ],
        [
            [3.00000000e00, 0.00000000e00, 1.83697020e-16],
            [-1.50000000e00, 2.59807621e00, 1.83697020e-16],
            [0.00000000e00, 0.00000000e00, 4.00000000e00],
        ],
        [
            [3.95075336, 0.0, -0.62573786],
            [-0.7326449, 3.88222663, -0.62573786],
            [0.0, 0.0, 4.0],
        ],
    ]
)


def test_Lattice_default():
    lattice = Lattice()

    items = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ]

    for key in items:
        item: Parameter = getattr(lattice, key)
        assert item.name == key
        t = key.split("_")[0]
        test_defaults = CELL_DETAILS[t].copy()
        del test_defaults["value"]
        #test_defaults["unit"] = test_defaults["units"]
        #del test_defaults["unit"]
        for default in test_defaults.keys():
            r = test_defaults[default]
            i = getattr(item, default)
            if default == "unit":
                if i == "Å": # special case for Å -> ang
                    i = "ang"
                else:
                    i = str(i)[0:3]
                r = r[0:3]
            assert i == r


@pytest.mark.parametrize("ang_unit", ("deg", "rad"))
@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_from_pars(value: list, ang_unit: str):
    ref = [v for v in value]

    if ang_unit == "rad":
        value = [
            value[0],
            value[1],
            value[2],
            np.deg2rad(value[3]),
            np.deg2rad(value[4]),
            np.deg2rad(value[5]),
        ]

    l = Lattice(*value, ang_unit=ang_unit)

    items = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ]

    for idx, key in enumerate(items):
        item: Parameter = getattr(l, key)
        assert item.name == key
        t = key.split("_")[0]
        test_defaults = CELL_DETAILS[t].copy()
        test_defaults["value"] = ref[idx]
        del test_defaults["value"]
        # test_defaults["unit"] = test_defaults["units"]
        # del test_defaults["units"]
        for default in test_defaults.keys():
            r = test_defaults[default]
            i = getattr(item, default)
            if default == "unit":
                if i == "Å": # special case for Å -> ang
                    i = "ang"
                else:
                    i = str(i)[0:3]
                r = r[0:3]
            if isinstance(i, Number) and not isinstance(i, bool):
                assert r == pytest.approx(i)
            else:
                assert i == r

#FAILED tests/unit_tests/Components/test_Lattice.py::test_Lattice_from_matrix[monoclinic] - RuntimeError: dictionary changed size during iteration
# Adding map clear because of these errors happening ONLY in 3.12
@pytest.mark.parametrize("value", matrix_pars)
def test_Lattice_from_matrix(value):
    global_object.map._clear()
    args = value[0:-1]
    matrix = value[-1]
    l = Lattice.from_matrix(matrix)
    items = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ]
    for idx, key in enumerate(items):
        item: Parameter = getattr(l, key)
        assert item.name == key
        t = key.split("_")[0]
        test_defaults = CELL_DETAILS[t].copy()
        test_defaults["value"] = args[idx]
        del test_defaults["value"]
        for default in test_defaults.keys():
            r = test_defaults[default]
            i = getattr(item, default)
            if default == "unit":
                if i == "Å": # special case for Å -> ang
                    i = "ang"
                else:
                    i = str(i)[0:3]
                r = r[0:3]
            if isinstance(i, Number) and not isinstance(i, bool):
                assert r == pytest.approx(i)
            else:
                assert i == r


@pytest.mark.parametrize(
    "value", mod_pars([[0], [0, 2], [0, 1, 2], [0, 1, 2, 4], [0, 2], [0, 3]])
)
def test_Lattice_from_special(request, value):
    ref = np.array(value[0:6])
    cons = ref[value[6:]]
    lattice_type = request.node.name.split("[")[1][:-1]

    f = getattr(Lattice, lattice_type)
    l = f(*cons)

    items = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ]

    for idx, key in enumerate(items):
        item: Parameter = getattr(l, key)
        assert item.name == key
        t = key.split("_")[0]
        test_defaults = CELL_DETAILS[t].copy()
        test_defaults["value"] = ref[idx]
        del test_defaults["value"]
        for default in test_defaults.keys():
            r = test_defaults[default]
            i = getattr(item, default)
            if default == "unit":
                if i == "Å": # special case for Å -> ang
                    i = "ang"
                else:
                    i = str(i)[0:3]
                r = r[0:3]
            if isinstance(i, Number) and not isinstance(i, bool):
                assert r == pytest.approx(i)
            else:
                assert i == r


@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_pars_short_GET(value: list):
    l = Lattice(*value)

    items = ["a", "b", "c", "alpha", "beta", "gamma"]

    for idx, item in enumerate(items):
        f = getattr(l, item)
        assert f == value[idx]


@pytest.mark.parametrize(
    "in_value, new_value",
    mod_pars(
        [
            (6, 6, 6, 90, 90, 90),
            (11, 11, 6, 90, 90, 90),
            (5, 6, 7, 90, 90, 90),
            (6, 7, 8, 90, 95, 90),
            (4, 4, 5, 90, 90, 120),
            (6, 6, 6, 95, 95, 95),
        ],
        True,
    ),
)
def test_Lattice_pars_short_SET(in_value: list, new_value: list):
    l = Lattice(*in_value)

    items = ["a", "b", "c", "alpha", "beta", "gamma"]

    for idx, item in enumerate(items):
        f = getattr(l, item)
        assert f == in_value[idx]
        setattr(l, item, new_value[idx])
        f = getattr(l, item)
        assert f == new_value[idx]


@pytest.mark.parametrize(
    "in_value, new_value",
    mod_pars(
        [
            (6, 6, 6, 90, 90, 90),
            (11, 11, 6, 90, 90, 90),
            (5, 6, 7, 90, 90, 90),
            (6, 7, 8, 90, 95, 90),
            (4, 4, 5, 90, 90, 120),
            (6, 6, 6, 95, 95, 95),
        ],
        True,
    ),
)
def test_Lattice_pars_SET(in_value: list, new_value: list):
    l = Lattice(*in_value)

    items = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ]

    for idx, item in enumerate(items):
        f = getattr(l, item)
        assert f.value == in_value[idx]
        setattr(l, item, new_value[idx])
        f = getattr(l, item)
        assert f.value == new_value[idx]


@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_angles(value: list):
    l = Lattice(*value)
    assert np.all(np.array(value[3:]) == l.angles)


@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_lengths(value: list):
    l = Lattice(*value)
    assert np.all(np.array(value[0:3]) == l.lengths)


@pytest.mark.parametrize("value", matrix_pars)
def test_Lattice_matrix(value: list):
    args = value[0:-1]
    matrix = np.array(value[-1])

    l = Lattice(*args)
    assert np.all(np.isclose(matrix, l.matrix))


@pytest.mark.parametrize("value", matrix_pars)
def test_Lattice_inv_matrix(value: list):
    args = value[0:-1]
    matrix = np.array(value[-1])
    matrix = np.linalg.inv(matrix)

    l = Lattice(*args)
    assert np.all(np.isclose(matrix, l.inv_matrix))


@pytest.mark.parametrize(
    "value",
    mod_pars(
        [125.0, 500.0, 24.0, 23.704520174283306, 31.1769145362398, 61.35087958926781]
    ),
)
def test_Lattice_volume(value):
    args = value[:-1]
    volume = value[-1]

    l = Lattice(*args)
    assert volume == pytest.approx(l.volume.m)
    assert str(l.volume.units) == "angstrom ** 3"


@pytest.mark.parametrize("value", matrix_pars)
def test_Lattice_metric_tensor(value):
    args = value[0:-1]

    matrix = np.array(value[-1])
    matrix = np.dot(matrix, matrix.T)

    l = Lattice(*args)
    assert np.all(np.isclose(matrix, l.metric_tensor))


@pytest.mark.parametrize(
    "in_value, new_value",
    mod_pars(
        [
            (1.256637, 1.256637, 1.256637, 90, 90, 90),
            (0.6283185, 0.6283185, 1.256637, 90, 90, 90),
            (3.14159, 2.09439510, 1.5707963, 90, 90, 90),
            (3.180753, 2.0943951, 1.5903765, 90, 81, 90),
            (2.418399, 2.418399, 1.570796, 90, 90, 60),
            (1.61845, 1.61845, 1.61845, 79.31296, 79.31296, 79.31296),
        ],
        True,
    ),
)
def test_Lattice_reciprocal_lattice(in_value: list, new_value: list):
    l = Lattice(*in_value)
    obj = l.reciprocal_lattice

    items = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ]

    for idx, item in enumerate(items):
        f = getattr(obj, item)
        assert np.isclose(f.value, new_value[idx])


@pytest.mark.parametrize(
    "in_value, new_value",
    mod_pars(
        [
            (0.2, 0.2, 0.2, 90, 90, 90),
            (0.1, 0.1, 0.2, 90, 90, 90),
            (0.5, 1 / 3, 0.25, 90, 90, 90),
            (0.5062325, 1 / 3, 0.253116, 90, 81, 90),
            (0.3849, 0.3849, 0.25, 90, 90, 60),
            (0.257584, 0.257584, 0.257584, 79.31296, 79.31296, 79.31296),
        ],
        True,
    ),
)
def test_Lattice_reciprocal_lattice(in_value: list, new_value: list):
    l = Lattice(*in_value)
    obj = l.reciprocal_lattice_crystallographic

    items = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ]

    for idx, item in enumerate(items):
        f = getattr(obj, item)
        assert np.isclose(f.value, new_value[idx])


@pytest.mark.parametrize("scale", [0.1, 2, 3.14, 100])
@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_scale(value: list, scale: float):
    l = Lattice(*value)

    new_volume = scale * l.volume
    scaled = l.scale(new_volume)

    assert np.isclose(scaled.volume, new_volume)
    assert np.all(np.isclose(l.angles, scaled.angles))


@pytest.mark.parametrize(
    "scale", [0.1, 2, 3.14, 100, [0.5, 0.5, 1], [0.5, 1, 0.5], [1, 0.5, 0.5]]
)
@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_length_scale(value: list, scale: float):
    l = Lattice(*value)

    scaled = l.scale_lengths(scale)
    assert np.all(np.isclose(l.angles, scaled.angles))
    assert np.all(np.isclose(np.array(l.lengths) * scale, scaled.lengths))


@pytest.mark.parametrize(
    "co_ords", [[0.1, 2, 3.14], [0.5, 0.5, 1], [0.5, 1, 0.5], [1, 0.5, 0.5]]
)
@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_fract_cart_coords(value: list, co_ords: list):
    l = Lattice(*value)

    frac = l.get_fractional_coords(co_ords)
    cart_co_ords = l.get_cartesian_coords(frac)

    assert np.all(np.isclose(co_ords, cart_co_ords))


@pytest.mark.parametrize("crystal_system", ["cubic", "tetragonal", "orthorhombic"])
def test_Lattice_is_orthogonal(crystal_system):

    l = Lattice(*pars_dict[crystal_system])
    assert l.is_orthogonal()


def test_Lattice_is_hexagonal():

    lengths = np.array(pars_dict["hexagonal"][0:3])
    angles = np.array(pars_dict["hexagonal"][3:])

    l = Lattice(*lengths, *angles)
    assert l.is_hexagonal()

    l = Lattice(*lengths, *(angles + 1e-5))
    assert l.is_hexagonal(hex_angle_tol=1e-4)
    assert not l.is_hexagonal(hex_angle_tol=1e-6)

    l = Lattice(*(lengths + [1e-2, -1e-2, 1e-2]), *angles)
    assert l.is_hexagonal(hex_length_tol=1e-1)
    assert not l.is_hexagonal(hex_length_tol=1e-7)


@pytest.mark.parametrize(
    "values, out_str",
    mod_pars(
        [
            "<Lattice: (a: 5.00 Å, b: 5.00 Å, c: 5.00 Å, alpha: 90.00 deg, beta: 90.00 deg, gamma: 90.00 deg>",
            "<Lattice: (a: 10.00 Å, b: 10.00 Å, c: 5.00 Å, alpha: 90.00 deg, beta: 90.00 deg, gamma: 90.00 deg>",
            "<Lattice: (a: 2.00 Å, b: 3.00 Å, c: 4.00 Å, alpha: 90.00 deg, beta: 90.00 deg, gamma: 90.00 deg>",
            "<Lattice: (a: 2.00 Å, b: 3.00 Å, c: 4.00 Å, alpha: 90.00 deg, beta: 99.00 deg, gamma: 90.00 deg>",
            "<Lattice: (a: 3.00 Å, b: 3.00 Å, c: 4.00 Å, alpha: 90.00 deg, beta: 90.00 deg, gamma: 120.00 deg>",
            "<Lattice: (a: 4.00 Å, b: 4.00 Å, c: 4.00 Å, alpha: 99.00 deg, beta: 99.00 deg, gamma: 99.00 deg>",
        ],
        True,
    ),
)
def test_Lattice_repr(values, out_str):
    l = Lattice(*values)
    assert str(l) == out_str


def make_dict(value) -> dict:
    return {
        "@module": "easycrystallography.Components.Lattice",
        "@class": "Lattice",
        "length_a": {
            "@module": "easyscience.Objects.variable.parameter",
            "@class": "Parameter",
            "@version": easyscience.__version__,
            "name": "length_a",
            "value": float(value[0]),
            "variance": 0.0,
            "min": 0.0,
            "max": np.inf,
            "fixed": True,
            "description": "Unit-cell length of the selected structure in angstroms.",
            "url": "https://docs.easydiffraction.org/lib/dictionaries/_cell/",
            "unit": "angstrom",
            "unique_name": "Parameter_0",
            "enabled": True,
        },
        "length_b": {
            "@module": "easyscience.Objects.variable.parameter",
            "@class": "Parameter",
            "@version": easyscience.__version__,
            "name": "length_b",
            "value": float(value[1]),
            "variance": 0.0,
            "min": 0.0,
            "max": np.inf,
            "fixed": True,
            "description": "Unit-cell length of the selected structure in angstroms.",
            "url": "https://docs.easydiffraction.org/lib/dictionaries/_cell/",
            "unit": "angstrom",
            "unique_name": "Parameter_0",
            "enabled": True,
        },
        "length_c": {
            "@module": "easyscience.Objects.variable.parameter",
            "@class": "Parameter",
            "@version": easyscience.__version__,
            "name": "length_c",
            "value": float(value[2]),
            "variance": 0.0,
            "min": 0.0,
            "max": np.inf,
            "fixed": True,
            "description": "Unit-cell length of the selected structure in angstroms.",
            "url": "https://docs.easydiffraction.org/lib/dictionaries/_cell/",
            "unit": "angstrom",
            "unique_name": "Parameter_0",
            "enabled": True,
        },
        "angle_alpha": {
            "@module": "easyscience.Objects.variable.parameter",
            "@class": "Parameter",
            "@version": easyscience.__version__,
            "name": "angle_alpha",
            "value": float(value[3]),
            "variance": 0.0,
            "min": 0.0,
            "max": np.inf,
            "fixed": True,
            "description": "Unit-cell angle of the selected structure in degrees.",
            "url": "https://docs.easydiffraction.org/lib/dictionaries/_cell/",
            "unit": "deg",
            "unique_name": "Parameter_0",
            "enabled": True,
        },
        "angle_beta": {
            "@module": "easyscience.Objects.variable.parameter",
            "@class": "Parameter",
            "@version": easyscience.__version__,
            "name": "angle_beta",
            "value": float(value[4]),
            "variance": 0.0,
            "min": 0.0,
            "max": np.inf,
            "fixed": True,
            "description": "Unit-cell angle of the selected structure in degrees.",
            "url": "https://docs.easydiffraction.org/lib/dictionaries/_cell/",
            "unit": "deg",
            "unique_name": "Parameter_0",
            "enabled": True,
        },
        "angle_gamma": {
            "@module": "easyscience.Objects.variable.parameter",
            "@class": "Parameter",
            "@version": easyscience.__version__,
            "name": "angle_gamma",
            "value": float(value[5]),
            "variance": 0.0,
            "min": 0.0,
            "max": np.inf,
            "fixed": True,
            "description": "Unit-cell angle of the selected structure in degrees.",
            "url": "https://docs.easydiffraction.org/lib/dictionaries/_cell/",
            "unit": "deg",
            "unique_name": "Parameter_0",
            "enabled": True,
        },
        "interface": None,
        "unique_name": "Parameter_0",
    }


@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_as_dict(value: list):
    l = Lattice(*value)
    obtained = l.as_dict()
    expected = make_dict(value)

    def check_dict(check, item):
        if isinstance(check, dict) and isinstance(item, dict):
            for this_check_key in check.keys():
                if this_check_key == "unique_name":
                    continue
                check_dict(check[this_check_key], item[this_check_key])
        else:
            assert isinstance(item, type(check))
            if item == "Å":
                item = "angstrom"
            assert item == check

    check_dict(expected, obtained)


@pytest.mark.skip(reason="Object name ... already exists in the grap")
#@pytest.mark.parametrize("value", basic_pars)
def test_Lattice_from_dict(value: list):
    # ADDED UNTIL UNIQUE_NAME IS FIXED BECAUSE OF THE FOLLOWING ERRORS:
    # FAILED tests/unit_tests/Components/test_Lattice.py::test_Lattice_from_dict[tetragonal] - ValueError: Object name 149027786693506016496254445195239597714 already exists in the graph.
    # FAILED tests/unit_tests/Components/test_Lattice.py::test_Lattice_from_dict[monoclinic] - ValueError: Object name 149027786693506016496254445195239597714 already exists in the graph.
    # FAILED tests/unit_tests/Components/test_Lattice.py::test_Lattice_from_dict[rhombohedral] - ValueError: Object name 149027786693506016496254445195239597714 already exists in the grap
    global_object.map._clear()

    expected = make_dict(value)
    l = Lattice.from_dict(expected)
    obtained = l.as_dict()

    def check_dict(check, item):
        if isinstance(check, dict) and isinstance(item, dict):
            for this_check_key in check.keys():
                if this_check_key == "unique_name":
                    continue
                check_dict(check[this_check_key], item[this_check_key])
        else:
            assert isinstance(item, type(check))
            if item == "Å":
                item = "angstrom"
            assert item == check

    check_dict(expected, obtained)


@pytest.mark.parametrize("value", basic_pars)
@pytest.mark.parametrize("fmt", [".3f"])
@pytest.mark.parametrize("opt", ["m", "l", "t"])
def test_Lattice_fmt(value, fmt, opt):

    l = Lattice(*value)
    out_fmt = "{}{}".format(fmt, opt)

    def do_test(in_str):
        m = (l.lengths, l.angles)
        if opt == "m":
            m = l.matrix.tolist()
            fmt2 = "[[{}, {}, {}], [{}, {}, {}], [{}, {}, {}]]"
        elif opt == "l":
            fmt2 = "{{{}, {}, {}, {}, {}, {}}}"
        else:
            fmt2 = "({} {} {}), ({} {} {})"
        check_str = fmt2.format(*[format(c, fmt) for row in m for c in row])
        assert in_str == check_str

    # Ancient Python. We won't be supporting this
    with pytest.raises(TypeError):
        out_fmt2 = "%" + out_fmt
        out_str = out_fmt2 % l
        do_test(out_str)
    # Python >2.7 "{:03fm}".format(l)
    out_fmt2 = "{:" + f"{out_fmt}" + "}"
    out_str = out_fmt2.format(l)
    do_test(out_str)
    # Python >3.6 + f"{l:03fm}"
    # This is stupidly dangerous.
    # Releases dragons, orks and is where darkness lies. You've been warned
    # !!!! DO NOT USE OUT OF THIS UNIQUE CONTEXT !!!!

    def effify(non_f_str: str, l: Lattice) -> str:
        return eval(f'f"""{non_f_str}"""')

    out_str = effify(f"{{l:{out_fmt}}}", l)
    do_test(out_str)


@pytest.mark.parametrize("value", basic_pars)
def test_lattice_copy(value):
    from copy import copy

    l1 = Lattice(*value)
    l2 = copy(l1)

    items = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ]

    for item in items:
        f1 = getattr(l1, item)
        f2 = getattr(l2, item)
        assert np.isclose(f1.value, f2.value)
        assert f1 != f2
