# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

from __future__ import annotations

from typing import TYPE_CHECKING
from typing import ClassVar
from typing import List
from typing import NoReturn
from typing import Tuple

from easycrystallography.Components.SpaceGroup import SpaceGroup as _SpaceGroup
from easycrystallography.Symmetry.SymOp import SymmOp

from .template import CIF_Template
from .template import gemmi

if TYPE_CHECKING:
    from easyscience.Utils.typing import B


class SpaceGroup(CIF_Template):
    _CIF_SECTION_NAME: ClassVar[str] = '_space_group'
    _CIF_CONVERSIONS: ClassVar[List[Tuple[str, str]]] = [
        ('space_group_HM_name', '_name_H-M_ref'),
        ('setting', '_IT_coordinate_system_code'),
        ('symmetry_ops', "'_symop.operation_xyz'"),
    ]
    _CIF_ALTERNATES: ClassVar[List[Tuple[str, List[str]]]] = [
        [
            'space_group_HM_name',
            [
                '_name_H-M_full',
                '_IT_number',
                '_name_Hall',
                '_name_H-M_alt',
                '.name_H-M_full',
                '.IT_number',
                '.name_Hall',
                '.name_H-M_alt',
            ],
        ],
        ['setting', ['_IT_coordinate_system_code', '.IT_coordinate_system_code']],
    ]

    def __init__(self, reference_class=_SpaceGroup):
        super().__init__()
        self._CIF_CLASS = reference_class

    def from_cif_block(self, block: gemmi.cif.Block) -> B:
        kwargs = {}
        for item in self._CIF_CONVERSIONS[0:2]:
            value = block.find_pair_item(self._CIF_SECTION_NAME + item[1])
            if value is None:
                a, b = zip(*self._CIF_ALTERNATES)
                if item[0] in a:
                    idx = a.index(item[0])
                    for i in b[idx]:
                        value = block.find_pair_item(self._CIF_SECTION_NAME + i)
                        if value is not None:
                            break
            if value is None:
                continue
            V, E, F = self.string_to_variable(value.pair[1])
            kwargs[item[0]] = V
        if not kwargs:
            item = self._CIF_CONVERSIONS[2]
            loop = list(block.find_values(self._CIF_SECTION_NAME + item[1]))
            ops = []
            for this_item in list(loop):
                ops.append(SymmOp.from_xyz_string(this_item))
            kwargs[item[0]] = ops
        return self._CIF_CLASS(**kwargs)

    def add_to_cif_block(self, obj: B, block: gemmi.cif.Block) -> NoReturn:
        if not obj.is_custom:
            for item in self._CIF_CONVERSIONS[0:2]:
                value = getattr(obj, item[0])
                if value:
                    block.set_pair(self._CIF_SECTION_NAME + item[1], self.variable_to_string(value))
        else:
            loop = block.init_loop('_space_group_symop.', ['id', 'operation_xyz'])
            for i, op in enumerate(obj.symmetry_ops):
                loop.add_row([str(i + 1), "'" + op.as_xyz_string().replace(' ', '') + "'"])

    def from_cif_string(self, cif_string: str) -> List[B]:
        if 'data_' not in cif_string:
            cif_string = 'data_temp\n' + cif_string

        cif_blocks = gemmi.cif.read_string(cif_string)
        objs = []
        for block in cif_blocks:
            objs.append(self.from_cif_block(block))
        return objs
