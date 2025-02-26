# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

from __future__ import annotations

from typing import TYPE_CHECKING
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional
from typing import Union

import numpy as np
from easyscience.Objects.Groups import BaseCollection
from easyscience.Objects.ObjectClasses import BaseObj
from easyscience.Objects.variable import DescriptorStr as Descriptor
from easyscience.Objects.variable import Parameter

from easycrystallography.Components.Lattice import Lattice
from easycrystallography.Components.Lattice import PeriodicLattice
from easycrystallography.Components.Site import Atoms
from easycrystallography.Components.Site import PeriodicAtoms
from easycrystallography.Components.Site import Site
from easycrystallography.Components.SpaceGroup import SpaceGroup
from easycrystallography.io.parser import Parsers

if TYPE_CHECKING:
    from easyscience.Utils.typing import iF


class Phase(BaseObj):
    _SITE_CLASS = Site
    _ATOMS_CLASS = Atoms

    cell = ClassVar[PeriodicLattice]
    _spacegroup = ClassVar[SpaceGroup]
    atoms = ClassVar[PeriodicAtoms]
    scale = ClassVar[Parameter]

    _REDIRECT = {
        'space_group': lambda obj: getattr(obj, '_spacegroup'),
    }

    def __init__(
        self,
        name: str,
        space_group: Optional[Union[SpaceGroup, str]] = None,
        cell: Optional[Union[Lattice, PeriodicLattice]] = None,
        atoms: Optional[Atoms] = None,
        scale: Optional[Parameter] = None,
        interface: Optional[iF] = None,
        enforce_sym: bool = True,
    ):
        self.name = name
        if space_group is None:
            space_group = SpaceGroup()
        if cell is None:
            cell = Lattice()
        if isinstance(cell, Lattice):
            cell = PeriodicLattice.from_lattice_and_spacegroup(cell, space_group)
        if atoms is None:
            atoms = Atoms('atoms')
        if scale is None:
            scale = Parameter('scale', 1, min=0)

        super(Phase, self).__init__(name, cell=cell, _spacegroup=space_group, atoms=atoms, scale=scale)
        if not enforce_sym:
            self.cell.clear_sym()
        self._enforce_sym = enforce_sym
        self.interface = interface

        self._extent = np.array([1, 1, 1])
        self._centre = np.array([0, 0, 0])
        self.atom_tolerance = 1e-4

    def add_atom(self, *args, **kwargs):
        """
        Add an atom to the crystal
        """
        supplied_atom = False
        for arg in args:
            if issubclass(arg.__class__, self._SITE_CLASS):
                self.atoms.append(arg)
                supplied_atom = True
        if not supplied_atom:
            atom = Site(*args, **kwargs)
            self.atoms.append(atom)

    def remove_atom(self, key):
        del self.atoms[key]

    def all_orbits(self, extent=None, magnetic_only: bool = False) -> Dict[str, np.ndarray]:
        """
        Generate all atomic positions from the atom array and symmetry operations over an extent.

        :return:  dictionary with keys of atom labels, containing numpy arrays of unique points in the extent
        (0, 0, 0) -> obj.extent
        :rtype: Dict[str, np.ndarray]
        """

        if extent is None:
            extent = self._extent

        offsets = np.array(
            np.meshgrid(
                range(0, extent[0] + 1),
                range(0, extent[1] + 1),
                range(0, extent[2] + 1),
            )
        ).T.reshape(-1, 3)

        orbits = self.get_orbits(magnetic_only=magnetic_only)
        for orbit_key in orbits.keys():
            orbit = orbits[orbit_key]
            site_positions = np.apply_along_axis(np.add, 1, offsets, orbit).reshape((-1, 3)) - self.center
            orbits[orbit_key] = (
                site_positions[
                    np.all(site_positions >= -self.atom_tolerance, axis=1)
                    & np.all(site_positions <= extent + self.atom_tolerance, axis=1),
                    :,
                ]
                + self.center
            )
        return orbits

    def get_orbits(self, magnetic_only: bool = False) -> Dict[str, np.ndarray]:
        """
        Generate all atomic positions from the atom array and symmetry operations over an extent.

        :return:  dictionary with keys of atom labels, containing numpy arrays of unique points in the extent
        (0, 0, 0) -> obj.extent
        :rtype: Dict[str, np.ndarray]
        """
        atoms = PeriodicAtoms.from_atoms(self.cell, self.atoms)
        orbits = atoms.get_orbits(magnetic_only=magnetic_only)
        return orbits

    @property
    def enforce_sym(self):
        return self._enforce_sym

    @enforce_sym.setter
    def enforce_sym(self, value: bool):
        if value:
            self.cell.enforce_sym()
        else:
            self.cell.clear_sym()

    @property
    def space_group(self):
        return self._spacegroup

    @space_group.setter
    def space_group(self, value: Union[SpaceGroup, str]):
        if isinstance(value, str):
            value = SpaceGroup(value)
        self._spacegroup = value

    def set_spacegroup(self, value):
        if self._enforce_sym:
            self.cell.space_group_HM_name = value
        else:
            self._spacegroup.space_group_HM_name = value

    @property
    def extent(self) -> np.ndarray:
        """
        Get the current extent in unit cells

        :return: current extent in unit cells
        :rtype: np.ndarray
        """
        return self._extent

    @extent.setter
    def extent(self, new_extent: Union[list, np.ndarray]):
        """
        The current extent of in unit cells. Default (1, 1, 1)

        :param new_extent: The new extent in unit cells.
        :type new_extent: Union[list, tuple, np.ndarray]
        """
        if isinstance(new_extent, list):
            new_extent = np.array(new_extent)
        if np.prod(new_extent.shape) != 3:
            raise ValueError
        new_extent = new_extent.reshape((3,))
        self._extent = new_extent

    @property
    def center(self) -> np.ndarray:
        """
        Get the center position

        :return: center position
        :rtype: np.ndarray
        """
        return self._centre

    @center.setter
    def center(self, new_center: Union[list, tuple, np.ndarray]):
        """
        Set the center position. Default (0, 0, 0)

        :param new_center: New center position.
        :type new_center: Union[list, tuple, np.ndarray]
        """

        if isinstance(new_center, list):
            new_center = np.array(new_center)
        if np.prod(new_center.shape) != 3:
            raise ValueError
        new_center = new_center.reshape((3,))
        self._centre = new_center

    def _generate_positions(self, site, extent) -> np.ndarray:
        """
        Generate all orbits for a given fractional position.
        """
        sym_op = self.space_group._sg_data.get_orbit
        offsets = np.array(
            np.meshgrid(
                range(0, extent[0] + 1),
                range(0, extent[1] + 1),
                range(0, extent[2] + 1),
            )
        ).T.reshape(-1, 3)
        return np.apply_along_axis(np.add, 1, offsets, np.array(sym_op(site.fract_coords))).reshape((-1, 3))

    def all_sites(self, extent=None) -> Dict[str, np.ndarray]:
        """
        Generate all atomic positions from the atom array and symmetry operations over an extent.
        :return:  dictionary with keys of atom labels, containing numpy arrays of unique points in the extent
        (0, 0, 0) -> obj.extent
        :rtype: Dict[str, np.ndarray]
        """
        if self.space_group is None:
            return {atom.label: atom.fract_coords for atom in self.atoms}

        if extent is None:
            extent = self._extent

        sites = {}
        for site in self.atoms:
            unique_sites = self._generate_positions(site, extent)
            site_positions = unique_sites - self.center
            sites[site.label.value] = (
                site_positions[
                    np.all(site_positions >= -self.atom_tolerance, axis=1)
                    & np.all(site_positions <= extent + self.atom_tolerance, axis=1),
                    :,
                ]
                + self.center
            )
        return sites

    @property
    def cif(self) -> str:
        s = ''
        cif_str_parser = Parsers('cif_str')
        with cif_str_parser.writer() as r:
            s += r.structure(self)
        s = s.strip()  # remove leading and trailing newlines
        return s

    @classmethod
    def from_cif_file(cls, filename):
        s = None
        with Parsers('cif').reader(filename) as r:
            s = r.structures(phase_class=cls)
        return s

    @classmethod
    def from_cif_string(cls, cif_string: str):
        s = None
        with Parsers('cif_str').reader() as r:
            s = r.structures(cif_string, phase_class=cls)
        return s


class Phases(BaseCollection):
    _SITE_CLASS = Site
    _ATOM_CLASS = Atoms
    _PHASE_CLASS = Phase

    def __init__(self, name: str = 'phases', *args, interface: Optional[iF] = None, **kwargs):
        """
        Generate a collection of crystals.

        :param name: Name of the crystals collection
        :type name: str
        :param args: objects to create the crystal
        :type args: *Phase
        """
        if not isinstance(name, str):
            raise AttributeError('Name should be a string!')

        super(Phases, self).__init__(name, *args, **kwargs)
        self.interface = interface

    def __repr__(self) -> str:
        return f'Collection of {len(self)} phases: {self.phase_names}'

    def __getitem__(self, idx: Union[int, slice]) -> Union[Parameter, Descriptor, BaseObj, BaseCollection]:
        if isinstance(idx, str) and idx in self.phase_names:
            idx = self.phase_names.index(idx)
        return super(Phases, self).__getitem__(idx)

    def __delitem__(self, key):
        if isinstance(key, str) and key in self.phase_names:
            key = self.phase_names.index(key)
        return super(Phases, self).__delitem__(key)

    def append(self, item: Phase):
        if not isinstance(item, Phase):
            raise TypeError('Item must be a Phase')
        if item.name in self.phase_names:
            raise AttributeError(f'A phase of name {item.name} already exists.')
        super(Phases, self).append(item)

    @property
    def phase_names(self) -> List[str]:
        return [phase.name for phase in self]

    @property
    def cif(self) -> str:
        s = ''
        cif_str_parser = Parsers('cif_str')
        with cif_str_parser.writer() as r:
            s += r.structures(self)
        return s

    @classmethod
    def from_cif_file(cls, filename):
        s = None
        with Parsers('cif').reader(filename) as r:
            s = r.structures(phases_class=cls, phase_class=cls._PHASE_CLASS)
        return s

    @classmethod
    def from_cif_string(cls, cif_string: str):
        s = None
        with Parsers('cif_str').reader() as r:
            s = r.structures(cif_string, phase_class=cls._PHASE_CLASS)
        return s
