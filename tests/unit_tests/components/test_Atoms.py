# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

from typing import List

import pytest
import numpy as np
from easycrystallography.Components.Site import Atoms, Site, _SITE_DETAILS

site_details = [Site('Al', 'Al'), Site('Fe', 'Fe3+'), Site('TEST', 'H')]


def gen_sites() -> List[Site]:
    args = []
    for edx in range(len(site_details)):
        key = f'from_{edx+1}_sites'
        args.append(pytest.param((site_details[0:(edx + 1)]), id=key))
    return args


def test_Atoms_empty():
    name = 'test'
    atoms = Atoms(name)
    assert len(atoms) == 0
    assert atoms.name == name


def test_Atoms_wrong_name():
    name = False
    with pytest.raises(TypeError):
        atoms = Atoms(name)


# def test_Atoms_wrong_args():
#     name = 'test'
#     args = [False, 'testing', (), []]
#     with pytest.raises(TypeError):
#         atoms = Atoms(name, *args)

@pytest.mark.parametrize('sites', gen_sites())
def test_Atoms_creation_args(sites: List[Site]):

    name = 'test'
    atoms = Atoms(name, *sites)
    assert len(atoms) == len(sites)
    assert atoms.name == name


@pytest.mark.parametrize('sites', gen_sites())
def test_Atoms_creation_dict(sites: List[Site]):

    name = 'test'
    d = {}
    for idx, site in enumerate(sites):
        d[str(idx)] = site
    atoms = Atoms(name, **d)
    assert len(atoms) == len(sites)
    assert atoms.name == name


@pytest.mark.parametrize('sites', gen_sites())
def test_Atoms_repr(sites: List[Site]):
    name = 'test'
    atoms = Atoms(name, *sites)
    assert str(atoms) == f'Collection of {len(atoms)} sites.'


@pytest.mark.parametrize('sites', gen_sites())
def test_Atoms_get_item_int(sites: List[Site]):
    name = 'test'
    atoms = Atoms(name, *sites)

    for atom, site in zip(atoms, sites):
        assert atom.label == site.label

    for idx, site in enumerate(sites):
        atom = atoms[idx]
        assert atom.label == site.label


@pytest.mark.parametrize('sites', gen_sites())
def test_Atoms_get_item_str(sites: List[Site]):
    name = 'test'
    atoms = Atoms(name, *sites)

    for site in sites:
        atom = atoms[site.label.value]
        assert atom.label == site.label