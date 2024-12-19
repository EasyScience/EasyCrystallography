# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING
from typing import ClassVar
from typing import List
from typing import NoReturn
from typing import Optional
from typing import Tuple
from typing import Union

import gemmi
import numpy as np
from easyscience.Objects.ObjectClasses import BaseObj
from easyscience.Objects.ObjectClasses import Descriptor

from easycrystallography.Symmetry.functions import get_default_it_coordinate_system_code_by_it_number
from easycrystallography.Symmetry.SymOp import SymmOp

SG_DETAILS = {
    'space_group_HM_name': {
        'name': 'hermann_mauguin',
        'description': 'Hermann-Mauguin symbols given in Table 4.3.2.1 of International Tables for Crystallography '
        'Vol. A (2002) or a Hermann-Mauguin symbol for a conventional or unconventional setting.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_space_group/',
        'value': 'P 1',
    },
    'setting': {
        'name': 'coordinate-code',
        'description': 'A qualifier taken from the enumeration list identifying which setting in International Tables '
        'for Crystallography Volume A (2002) (IT) is used.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_space_group/',
        'value': '\x00',
    },
    'symmetry_ops': {
        'name': 'symmetry-ops',
        'description': 'A list of symmetry operations, each given as a 4x4 matrix in the form of a `SymmOp` object.',
        'value': [SymmOp.from_xyz_string('x,y,z')],
    },
}

if TYPE_CHECKING:
    import numpy.typing as npt
    from easyscience.Utils.typing import iF

    from easycrystallography.Components.Site import S

    T = Union[S, npt.ArrayLike]


_D_REDIRECT = deepcopy(Descriptor._REDIRECT)
_D_REDIRECT['value'] = lambda obj: ';'.join([r.as_xyz_string() for r in obj.raw_value.tolist()])


class easyOp(Descriptor):
    _REDIRECT = _D_REDIRECT


class SpaceGroup(BaseObj):
    _space_group_HM_name: ClassVar[Descriptor]
    _setting: ClassVar[Descriptor]
    _symmetry_ops: ClassVar[Descriptor]

    _REDIRECT = {'symmetry_ops': lambda obj: None if obj._sg_data is not None else obj._symmetry_ops}

    def __init__(
        self,
        space_group_HM_name: Optional[Descriptor, str] = None,
        setting: Optional[Descriptor, str] = None,
        symmetry_ops: Optional[List[SymmOp]] = None,
        interface: Optional[iF] = None,
    ):
        """
        Generate a spacegroup object from it's Hermann-Mauguin symbol and setting. The setting can be a part of the
        Hermann-Mauguin symbol.

        :param space_group_HM_name: Hermann-Mauguin symbol
        :param setting: Optional setting for the space group
        :param interface: Interface to the calculator
        """

        super(SpaceGroup, self).__init__(
            'space_group',
            _space_group_HM_name=Descriptor(**SG_DETAILS['space_group_HM_name']),
            _setting=Descriptor(**SG_DETAILS['setting']),
            _symmetry_ops=easyOp(**SG_DETAILS['symmetry_ops']),
        )

        if space_group_HM_name:
            self._space_group_HM_name = space_group_HM_name
        if setting is not None:
            self._setting = setting

        kwargs = {
            'new_spacegroup': self._space_group_HM_name.raw_value,
            'new_setting': self._setting.raw_value,
            'operations_set': None,
            'set_internal': True,
        }
        if symmetry_ops is not None:
            kwargs['operations_set'] = symmetry_ops
        self.__on_change(**kwargs)
        self.interface = interface
        self._cell = None

    @classmethod
    def from_int_number(cls, int_number: int, hexagonal=True, interface: Optional[iF] = None):
        """
        Generate a spacegroup object from it's spacegroup number (1-231).
        :param int_number: spacegroup number
        :param hexagonal: Should a hexagonal setting be used?
        :param interface: Interface to the calculator
        """
        sg = gemmi.find_spacegroup_by_number(int(int_number))
        setting = None
        if int_number in [146, 148, 155, 160, 161, 166, 167]:
            if hexagonal:
                setting = 'H'
            else:
                setting = 'R'
        return cls(sg.hm, setting, interface=interface)

    @classmethod
    def from_symOps(cls, sym_ops: List[SymmOp], interface: Optional[iF] = None):
        """
        Create a space group from a list of easycrystallography symmetry operations.

        :param sym_ops: List of easycrystallography symmetry operations
        :param interface: Interface to the calculator
        """
        ops = [gemmi.Op(op.as_xyz_string()) for op in sym_ops]
        sg_data = gemmi.find_spacegroup_by_ops(gemmi.GroupOps(ops))
        return cls(sg_data.hm, sg_data.ext, interface=interface)

    @classmethod
    def from_symMatrices(
        cls,
        rotations: List[np.ndarray],
        translations: List[np.ndarray],
        interface: Optional[iF] = None,
    ):
        """
        Create a space group from a lists of rotations and translations. The number of rotations and translations must
        be the same.

        :param rotations: Array of rotation matrices [n*[3x3]]
        :param translations: Array of translations [n*[1x3]]
        :param interface: Interface to the calculator
        """
        ops = []
        for rot, tran in zip(rotations, translations):
            ops.append(SymmOp.from_rotation_and_translation(rot, tran))
        return cls.from_symOps(ops, interface=interface)

    @classmethod
    def from_generators(
        cls,
        rotations: List[np.ndarray],
        translations: List[np.ndarray],
        interface: Optional[iF] = None,
    ):
        """
        Create a space group from a lists of rotations and translations. Each translation is applied to each rotation.

        :param rotations: Array of rotation matrices [n*[3x3]]
        :param translations: Array of translations [m*[1x3]]
        :param interface: Interface to the calculator
        """
        rots = len(translations) * rotations
        trans = [np.array(d) for d in np.array(translations).repeat(len(rotations), axis=0).tolist()]
        return cls.from_symMatrices(rots, trans, interface=interface)

    @classmethod
    def from_gemmi_operations(
        cls,
        operations: Union[gemmi.GroupOps, List[gemmi.Op]],
        interface: Optional[iF] = None,
    ):
        """
        Generate a space group from a list of gemmi operations or a gemmi group of operations.

        :param operations: Operations which define a space group
        :param interface: Interface to the calculator
        """
        ops = []
        if isinstance(operations, gemmi.GroupOps):
            ops = operations
        else:
            for op in operations:
                if isinstance(op, gemmi.Op):
                    ops.append(op)
                else:
                    raise TypeError('Operations must be of type gemmi.Op')
            ops = gemmi.GroupOps(ops)
        sg_data = gemmi.find_spacegroup_by_ops(ops)
        return cls(sg_data.hm, sg_data.ext, interface=interface)

    @classmethod
    def from_xyz_string(cls, xyz_string: Union[str, List[str]], interface: Optional[iF] = None):
        """
        Create a space group from a string or list of strings of the form 'x,y,z'.

        :param xyz_string: String defining space group operators
        :param interface: Interface to the calculator
        """
        if isinstance(xyz_string, str):
            xyz_string = xyz_string.split(';')
        ops = []
        for xyz in xyz_string:
            ops.append(SymmOp.from_xyz_string(xyz))
        return cls.from_symOps(ops, interface)

    @staticmethod
    def find_settings_by_number(number):
        """
        Find the IT_coordinate_system_code by group's number.
        gemmi doesn't do it natively.
        """
        ext = []
        for item in gemmi.spacegroup_table():
            if item.number > number:
                break
            if item.number == number:
                st = ''
                # Cases where ext and qualifier are not empty
                if item.ext and item.ext != '\x00':
                    st += item.ext
                if item.qualifier:
                    st += item.qualifier
                # special cases of defaul settings, not explicitly defined in gemmi
                system = item.crystal_system_str()
                if system == 'orthorhombic' and not item.qualifier:
                    st += 'abc'
                elif (system == 'trigonal' or system == 'hexagonal') and not item.qualifier:
                    st += 'h'
                # failed, just assign "1" to triclinic/monoclinic/tetragonal
                if not st:
                    st = '1'
                ext.append(st)
        return ext

    def __on_change(
        self,
        new_spacegroup: Union[int, str],
        new_setting: Optional[str] = None,
        operations_set: Optional[List[SymmOp]] = None,
        set_internal: bool = True,
    ) -> Tuple[str, str, List[SymmOp]]:
        """
        Internal function to update the space group. This function is called when the space group is changed. It checks
        the form of the imputs, generates reference data, and updates the internal data (if requested).

        :param new_spacegroup: New space group number or name
        :param new_setting: New space group setting
        :param set_internal: Should internal objects be updated
        """
        setting = '\x00'
        if operations_set is None:
            if isinstance(new_spacegroup, str):
                if ':' in new_spacegroup:
                    new_spacegroup, new_setting = new_spacegroup.split(':')
                sg_data = gemmi.find_spacegroup_by_name(new_spacegroup)
                if sg_data is None:
                    try:
                        sg_data = gemmi.find_spacegroup_by_ops(gemmi.symops_from_hall(new_spacegroup))
                    except RuntimeError:
                        sg_data = None
            else:
                sg_data = gemmi.find_spacegroup_by_number(int(new_spacegroup))

            if sg_data is None:
                raise ValueError(f"Spacegroup '{new_spacegroup}' not found in database.")

            settings = self.find_settings_by_number(sg_data.number)
            hm_name = sg_data.hm
            reference = get_default_it_coordinate_system_code_by_it_number(sg_data.number)

            if new_setting is None or new_setting == '' or new_setting == '\x00':
                if reference is not None:
                    setting = reference
            else:
                try:
                    new_setting = int(new_setting)
                except ValueError:
                    pass
                new_setting = str(new_setting)
                if new_setting != reference:
                    # modify the space group with the new setting
                    new_sg_data = gemmi.find_spacegroup_by_name(sg_data.hm + ':' + new_setting)
                    if new_sg_data is None and new_setting in settings:
                        # this can be because the setting is the "default" setting which gemmi treats as a blank
                        new_sg_data = gemmi.find_spacegroup_by_name(sg_data.hm + ':' + reference)
                    if new_sg_data is None:
                        raise ValueError(f"Spacegroup '{new_spacegroup}:{new_setting}' not found in database.")
                    sg_data = new_sg_data
                    setting = get_default_it_coordinate_system_code_by_it_number(sg_data.number)
                else:
                    setting = new_setting
                if new_setting in settings:
                    setting = new_setting

            if operations_set is None:
                operations_set = sg_data.operations()
            operations = [
                SymmOp.from_rotation_and_translation(np.array(op.rot) / op.DEN, np.array(op.tran) / op.DEN)
                for op in operations_set
            ]
        else:
            sg_data = None
            hm_name = 'custom'
            if isinstance(operations_set, str):
                operations_set = [SymmOp.from_xyz_string(s) for s in operations_set.split(';')]
            operations = operations_set
        if set_internal:
            self._sg_data = sg_data
            self._space_group_HM_name.value = hm_name
            self._setting.value = setting
            self._symmetry_ops.value = operations
        return sg_data, setting, operations

    @property
    def is_custom(self) -> bool:
        return self._sg_data is None

    @property
    def setting(self) -> Optional[Descriptor]:
        """
        Space group setting. If the space group does not have a setting, this will be None.

        :return: Space group setting
        """
        setting_str = self._setting.raw_value
        if setting_str == '\x00':
            return None  # no setting
        return self._setting

    @setting.setter
    def setting(self, new_setting: str) -> NoReturn:
        """
        Set the space group setting.

        :param new_setting: Space group setting
        """
        _, setting, _ = self.__on_change(self._space_group_HM_name.raw_value, new_setting, set_internal=True)

    @property
    def it_coordinate_system_code(self) -> Optional[Descriptor]:
        """
        Space group setting. If the space group does not have a setting, this will be None.
        Equivalent to setting, defined to satisfy CIF Template

        :return: Space group setting
        """
        setting_str = self._setting.raw_value
        if setting_str == '\x00':
            return None  # no setting
        return self._setting

    @it_coordinate_system_code.setter
    def it_coordinate_system_code(self, new_setting: str) -> NoReturn:
        """
        Set the space group setting.
        Equivalent to setting, defined to satisfy CIF Template

        :param new_setting: Space group setting
        """
        _, setting, _ = self.__on_change(self._space_group_HM_name.raw_value, new_setting, set_internal=True)

    @property
    def setting_str(self) -> str:
        """
        Space group setting as a string. If the space group does not have a setting, this will be an empty string.

        :return: Space group setting as a string
        """
        if self.setting is None:
            return ''
        return self._setting.raw_value

    @property
    def space_group_HM_name(self) -> Descriptor:
        """
        Space group name as defined by a Hermann-Mauguin symbol

        :return: Space group name as easyCore Descriptor
        """
        return self._space_group_HM_name

    @space_group_HM_name.setter
    def space_group_HM_name(self, value: str) -> NoReturn:
        """
        Set the space group name as defined by a Hermann-Mauguin symbol

        :param value: Space group name as a string
        """
        self.__on_change(value, set_internal=True)

    @property
    def name_hm_alt(self) -> Descriptor:
        """
        Space group name as defined by a Hermann-Mauguin symbol
        Equivalent to space_group_HM_name, defined to satisfy CIF Template

        :return: Space group name as easyCore Descriptor
        """
        return self._space_group_HM_name

    @name_hm_alt.setter
    def name_hm_alt(self, value: str) -> NoReturn:
        """
        Set the space group name as defined by a Hermann-Mauguin symbol
        Equivalent to space_group_HM_name, defined to satisfy CIF Template

        :param value: Space group name as a string
        """
        self.__on_change(value, set_internal=True)

    @property
    def hermann_mauguin(self) -> str:
        """
        Space group name as defined by a Hermann-Mauguin symbol

        :return: Space group name as a string
        """
        return self._space_group_HM_name.raw_value

    @property
    def name_hall(self) -> str:
        """
        Hall symbol of the space group
        Equivalent to hall_symbol, defined to satisfy CIF Template

        :return: Hall symbol of the space group
        """
        hall = None
        if not self.is_custom:
            hall = self._sg_data.hall
        return hall

    @property
    def hall_symbol(self) -> str:
        """
        Hall symbol of the space group

        :return: Hall symbol of the space group
        """
        hall = None
        if not self.is_custom:
            hall = self._sg_data.hall
        return hall

    @property
    def int_number(self) -> int:
        """
        International number of the space group

        :return: International number of the space group
        """
        n = None
        if not self.is_custom:
            n = self._sg_data.number
        return n

    @int_number.setter
    def int_number(self, new_it_number: int) -> NoReturn:
        """
        Set the spacegroup by its international number

        :param new_it_number: International number of the new space group
        """
        self.__on_change(new_it_number, set_internal=True)

    @property
    def it_number(self) -> int:
        """
        International number of the space group
        Equivalent to int_number, defined to satisfy CIF Template

        :return: International number of the space group
        """
        n = None
        if not self.is_custom:
            n = self._sg_data.number
        return n

    @it_number.setter
    def it_number(self, new_it_number: int) -> NoReturn:
        """
        Set the spacegroup by its international number
        Equivalent to int_number, defined to satisfy CIF Template

        :param new_it_number: International number of the new space group
        """
        self.__on_change(new_it_number, set_internal=True)

    @property
    def crystal_system(self) -> str:
        """
        Which crystal system the space group belongs to

        :return: Crystal system of the space group
        """
        s = ''
        if not self.is_custom:
            s = self._sg_data.crystal_system_str()
        return s

    @property
    def symmetry_ops(self) -> List[SymmOp]:
        """
        List of symmetry operations of the space group

        :return: List of symmetry operations of the space group
        """
        return self._symmetry_ops.raw_value

    @symmetry_ops.setter
    def symmetry_ops(self, new_ops: List[SymmOp]) -> NoReturn:
        """
        Set the symmetry operations of the space group

        :param new_ops: List of new symmetry operations
        """
        self.__on_change(
            self._space_group_HM_name.raw_value,
            operations_set=new_ops,
            set_internal=True,
        )

    @property
    def symmetry_xyz(self) -> str:
        """
        Symmetry operations of the space group as a string

        :return: String of symmetry operations of the space group
        """
        return ';'.join([op.as_xyz_string() for op in self._symmetry_ops.raw_value])

    @property
    def is_reference_setting(self) -> bool:
        """
        Is the space group setting the reference setting?

        :return: Is the space group setting the reference setting?
        """
        r = True
        if not self.is_custom:
            r = self._sg_data.is_reference_setting()
        return r

    def symmetry_matrices(self) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """
        Get the rotational and translational matrices of the space group

        :return: Rotation and translation matrices
        """
        ops = self.symmetry_ops
        return [op.rotation_matrix.copy() for op in ops], [op.translation_vector.copy() for op in ops]

    def get_orbit(self, point: T, tol: float = 1e-5) -> np.ndarray:
        """
        Returns the orbit for a point.

        :param point: Point to get the orbit for
        :param tol: Tolerance for the orbit
        :return: Orbits of the point
        """
        orbit = []
        if not hasattr(point, '__iter__'):
            point = point.fract_coords
        point = np.array(point)
        for o in self.symmetry_ops:
            pp = np.array(o.operate(point))
            if not in_array_list(orbit, pp, tol=tol):
                orbit.append(pp)
        return np.array(orbit)

    def get_site_multiplicity(self, site: T, tol=1e-5) -> int:
        """
        Get the multiplicity of a given site

        :param site: Site to get the multiplicity of
        :param tol: Tolerance for the orbit
        :return: Multiplicity of the site
        """
        if not hasattr(site, '__iter__'):
            site = site.fract_coords
        site = np.array(site)
        multiplicity = 1
        for o in self.symmetry_ops:
            new_site = o.operate(site)
            if np.isclose(new_site, site, atol=tol).all():
                multiplicity += 1
        return multiplicity

    def __repr__(self) -> str:
        out_str = "<Spacegroup: system: '{:s}', number: {}, H-M: '{:s}'".format(
            self.crystal_system, self.int_number, self.hermann_mauguin
        )
        if self.setting_str:
            out_str = "{:s} setting: '{:s}'".format(out_str, self.setting_str)
        return out_str + '>'


def in_array_list(array_list, a, tol=1e-5) -> bool:
    """
    Extremely efficient nd-array comparison using numpy's broadcasting. This
    function checks if a particular array a, is present in a list of arrays.
    It works for arrays of any size, e.g., even matrix searches.

    Args:
        array_list ([array]): A list of arrays to compare to.
        a (array): The test array for comparison.
        tol (float): The tolerance. Defaults to 1e-5. If 0, an exact match is
            done.

    Returns:
        (bool)
    """
    if len(array_list) == 0:
        return False
    axes = tuple(range(1, a.ndim + 1))
    if not tol:
        return np.any(np.all(np.equal(array_list, a[None, :]), axes))
    return np.any(np.sum(np.abs(array_list - a[None, :]), axes) < tol)
