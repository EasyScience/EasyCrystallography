# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

from __future__ import annotations

from typing import TYPE_CHECKING
from typing import ClassVar
from typing import Dict
from typing import List
from typing import NoReturn
from typing import Tuple

from easycrystallography.Components.AtomicDisplacement import AtomicDisplacement as _AtomicDisplacement
from easycrystallography.Components.Site import Atoms as _Atoms
from easycrystallography.Components.Susceptibility import MagneticSusceptibility as _MagneticSusceptibility

from .template import CIF_Template
from .template import gemmi

if TYPE_CHECKING:
    from easyscience.Utils.typing import B


class AtomicDisplacement(CIF_Template):
    _CIF_SECTION_NAME: ClassVar[str] = '_atom_site'

    _CIF_ADP_ISO_CONVERSIONS = [
        ('label', '_label', '.label'),
        ('adp_type', '_adp_type', '.adp_type'),
        ('Biso', '_B_iso_or_equiv', '.B_iso_or_equiv'),
        ('Uiso', '_U_iso_or_equiv', '.U_iso_or_equiv'),
    ]
    _CIF_ADP_ANISO_CONVERSIONS = [
        ('label', '_aniso_label'),
        ('B_11', '_aniso_B_11'),
        ('B_12', '_aniso_B_12'),
        ('B_13', '_aniso_B_13'),
        ('B_22', '_aniso_B_22'),
        ('B_23', '_aniso_B_23'),
        ('B_33', '_aniso_B_33'),
        ('U_11', '_aniso_U_11'),
        ('U_12', '_aniso_U_12'),
        ('U_13', '_aniso_U_13'),
        ('U_22', '_aniso_U_22'),
        ('U_23', '_aniso_U_23'),
        ('U_33', '_aniso_U_33'),
    ]

    def __init__(self, reference_class=_AtomicDisplacement):
        super().__init__()
        self._CIF_CLASS = reference_class

    def from_cif_block(self, block: gemmi.cif.Block) -> Dict[str, B]:
        atom_dict = {}
        # ADP CHECKER
        keys = [
            self._CIF_SECTION_NAME + name[1] if 'label' in name[1] else '?' + self._CIF_SECTION_NAME + name[1]
            for name in self._CIF_ADP_ISO_CONVERSIONS
        ]
        table = block.find(keys)
        if table.loop is None:
            # this means it might be the dictionary CIF
            keys = [
                self._CIF_SECTION_NAME + name[2] if 'label' in name[2] else '?' + self._CIF_SECTION_NAME + name[2]
                for name in self._CIF_ADP_ISO_CONVERSIONS
            ]
            table = block.find(keys)
        for row in table:
            kwargs = {}
            errors = {}
            is_fixed = {}
            if not row.has(0):
                continue
            if row.has(1):
                kwargs['adp_type'] = row[1]
            else:
                if row.has(2):
                    kwargs['adp_type'] = 'Biso'
                elif row.has(3):
                    kwargs['adp_type'] = 'Uiso'
                else:
                    continue

            for i in range(2, 4):
                if row.has(i):
                    V, E, F = self.string_to_variable(row[i])
                    kwargs[kwargs['adp_type']] = V
                    if E:
                        errors[kwargs['adp_type']] = E
                    if F is not None and not F:
                        is_fixed[kwargs['adp_type']] = F
            obj = _AtomicDisplacement(**kwargs)
            for error in errors.keys():
                setattr(getattr(obj, error), 'error', errors[error])
            for atr in is_fixed.keys():
                setattr(getattr(obj, atr), 'fixed', is_fixed[atr])
            atom_dict[row[0]] = {'adp': obj}

        keys = [
            self._CIF_SECTION_NAME + name[1] if 'label' in name[1] else '?' + self._CIF_SECTION_NAME + name[1]
            for name in self._CIF_ADP_ANISO_CONVERSIONS
        ]
        table = block.find(keys)
        for row in table:
            kwargs = {}
            errors = {}
            is_fixed = {}
            if not row.has(0):
                continue
            if row.has(1):
                kwargs['adp_type'] = 'Bani'
                idx = 1
            elif row.has(7):
                idx = 7
                kwargs['adp_type'] = 'Uani'
            else:
                continue
            for i in range(idx, idx + 6):
                ll = self._CIF_ADP_ANISO_CONVERSIONS[i][0]
                if row.has(i):
                    V, E, F = self.string_to_variable(row[i])
                    kwargs[ll] = V
                    if E:
                        errors[ll] = E
                    if F is not None and not F:
                        is_fixed[ll] = F
            obj = _AtomicDisplacement(**kwargs)
            for error in errors.keys():
                setattr(getattr(obj, error), 'error', errors[error])
            for atr in is_fixed.keys():
                setattr(getattr(obj, atr), 'fixed', is_fixed[atr])
            atom_dict[row[0]] = {'adp': obj}
        return atom_dict

    def add_to_cif_block(self, obj: B, block: gemmi.cif.Block) -> NoReturn:
        # Add the additional anisotropic loops
        lines = []
        names = [self._CIF_ADP_ANISO_CONVERSIONS[0][1]]
        values_offset = 1
        # find the adp type first
        for atom in obj:
            if not hasattr(atom, 'adp'):
                continue
            if self.variable_to_string(atom.adp.adp_type) == 'Uani':
                values_offset = 7
            break
        for i in range(values_offset, values_offset + 6):
            names.append(self._CIF_ADP_ANISO_CONVERSIONS[i][1])

        for atom in obj:
            if not hasattr(atom, 'adp'):
                continue
            if not self.variable_to_string(atom.adp.adp_type) == 'Uani':  # non-anisotropic ADP on this atom
                continue
            line = [self.variable_to_string(atom.__getattribute__(self._CIF_ADP_ANISO_CONVERSIONS[0][0]))]
            # names.append(self._CIF_ADP_ANISO_CONVERSIONS[values_offset+i][1])
            for keys in self._CIF_ADP_ANISO_CONVERSIONS[values_offset : values_offset + 6]:
                key, cif_key = keys
                s = self.variable_to_string(atom.adp.__getattribute__(key))
                line.append(s)
            lines.append(line)

        # _, names = zip(*self._CIF_ADP_ANISO_CONVERSIONS[1])

        loop = block.init_loop(self._CIF_SECTION_NAME, names)
        for line in lines:
            loop.add_row(line)

    def iso_adp_block(self, obj: B, block: gemmi.cif.Block):
        # Add the isotropic objects
        labels = []
        objs = []

        if len(obj) == 0:
            return labels, objs
        adp0 = getattr(obj[0], 'adp', None)
        if adp0 is None:
            return labels, objs
        adp_type = self.variable_to_string(adp0.adp_type)
        labels = self._CIF_ADP_ISO_CONVERSIONS.copy()
        # should we have `Uani` in the line or only `Uiso`?
        if adp_type[0].lower() == 'u':
            del labels[2]
        else:
            del labels[3]
        ################################
        objs = [[getattr(_obj, 'adp')] for _obj in obj]
        del labels[0]
        return [[labels]] * len(objs), objs

    def from_cif_string(self, cif_string: str) -> List[Dict[str, B]]:
        if 'data_' not in cif_string:
            cif_string = 'data_temp\n' + cif_string
        cif_blocks = gemmi.cif.read_string(cif_string)
        objs = []
        for block in cif_blocks:
            objs.append(self.from_cif_block(block))
        return objs


class MagneticSusceptibility(CIF_Template):
    _CIF_SECTION_NAME: ClassVar[str] = '_atom_site'
    _CIF_MSP_ANISO_CONVERSIONS = [
        ('label', '_susceptibility_label'),
        ('msp_type', '_susceptibility_chi_type'),
        ('chi_11', '_susceptibility_chi_11'),
        ('chi_12', '_susceptibility_chi_12'),
        ('chi_13', '_susceptibility_chi_13'),
        ('chi_22', '_susceptibility_chi_22'),
        ('chi_23', '_susceptibility_chi_23'),
        ('chi_33', '_susceptibility_chi_33'),
    ]

    def __init__(self, reference_class=_MagneticSusceptibility):
        super().__init__()
        self._CIF_CLASS = reference_class

    def from_cif_block(self, block: gemmi.cif.Block) -> Dict[str, B]:
        atom_dict = {}
        keys = [
            self._CIF_SECTION_NAME + name[1] if 'label' in name[1] else '?' + self._CIF_SECTION_NAME + name[1]
            for name in self._CIF_MSP_ANISO_CONVERSIONS
        ]
        table = block.find(keys)
        for row in table:
            kwargs = {}
            errors = {}
            is_fixed = {}
            if not row.has(0):
                continue
            if row.has(1):
                kwargs['msp_type'] = row[1]
            else:
                continue
            idx = 2
            for i in range(idx, idx + 6):
                ll = self._CIF_MSP_ANISO_CONVERSIONS[i][0]
                if row.has(i):
                    V, E, F = self.string_to_variable(row[i])
                    kwargs[ll] = V
                    if E:
                        errors[ll] = E
                    if F is not None and not F:
                        is_fixed[ll] = F
            obj = _MagneticSusceptibility(**kwargs)
            for error in errors.keys():
                setattr(getattr(obj, error), 'error', errors[error])
            for atr in is_fixed.keys():
                setattr(getattr(obj, atr), 'fixed', is_fixed[atr])
            atom_dict[row[0]] = {'msp': obj}
        return atom_dict

    def add_to_cif_block(self, obj: B, block: gemmi.cif.Block) -> NoReturn:
        # Then add the additional loops
        lines = []
        for atom in obj:
            if not hasattr(atom, 'msp'):
                continue
            line = [self.variable_to_string(atom.__getattribute__(self._CIF_MSP_ANISO_CONVERSIONS[0][0]))]
            for keys in self._CIF_MSP_ANISO_CONVERSIONS[1:]:
                key, cif_key = keys
                s = self.variable_to_string(atom.msp.__getattribute__(key))
                line.append(s)
            lines.append(line)
        _, names = zip(*self._CIF_MSP_ANISO_CONVERSIONS)
        loop = block.init_loop(self._CIF_SECTION_NAME, names)
        for line in lines:
            loop.add_row(line)

    def from_cif_string(self, cif_string: str) -> List[B]:
        if 'data_' not in cif_string:
            cif_string = 'data_temp\n' + cif_string
        cif_blocks = gemmi.cif.read_string(cif_string)
        objs = []
        for block in cif_blocks:
            objs.append(self.from_cif_block(block))
        return objs


class Atoms(CIF_Template):
    _CIF_SECTION_NAME: ClassVar[str] = '_atom_site'
    _CIF_CONVERSIONS: ClassVar[List[Tuple[str, str]]] = [
        ('label', '_label', '.label'),
        ('specie', '_type_symbol', '.type_symbol'),
        ('fract_x', '_fract_x', '.fract_x'),
        ('fract_y', '_fract_y', '.fract_y'),
        ('fract_z', '_fract_z', '.fract_z'),
        ('occupancy', '_occupancy', '.occupancy'),
    ]

    def __init__(self, reference_class=_Atoms):
        super().__init__()
        self._CIF_CLASS = reference_class

    def _site_runner(self, block):
        keys = [
            self._CIF_SECTION_NAME + name[1] if 'occupancy' not in name[1] else '?' + self._CIF_SECTION_NAME + name[1]
            for name in self._CIF_CONVERSIONS
        ]
        table = (
            block.find(keys)
            if (table := block.find(keys)).loop is not None
            else block.find(
                [
                    self._CIF_SECTION_NAME + name[2] if 'occupancy' not in name[2] else '?' + self._CIF_SECTION_NAME + name[2]
                    for name in self._CIF_CONVERSIONS
                ]
            )
        )
        atom_dict = {}
        error_dict = {}
        fixed_dict = {}
        for row in table:
            kwargs = {}
            errors = {}
            is_fixed = {}
            for idx, item in enumerate(self._CIF_CONVERSIONS):
                ec_name, _, _ = item
                if row.has(idx):
                    V, E, F = self.string_to_variable(row[idx])
                    kwargs[ec_name] = V
                    if E:
                        errors[ec_name] = E
                    if F is not None and not F:
                        is_fixed[ec_name] = F
            atom_dict[kwargs['label']] = kwargs
            error_dict[kwargs['label']] = errors
            fixed_dict[kwargs['label']] = is_fixed

        # ADP CHECKER
        ADP_RUNNER = AtomicDisplacement()
        adp_dict = ADP_RUNNER.from_cif_block(block)
        for label, adp in adp_dict.items():
            if label in atom_dict:
                atom_dict[label].update(adp)

        # MSP Checker
        MSP_RUNNER = MagneticSusceptibility()
        msp_dict = MSP_RUNNER.from_cif_block(block)
        for label, msp in msp_dict.items():
            if label in atom_dict:
                atom_dict[label].update(msp)

        atoms = []
        for a in atom_dict.values():
            obj = self._CIF_CLASS._SITE_CLASS(**a)
            label = obj.label.value
            if label in error_dict.keys():
                for atr in error_dict[label].keys():
                    setattr(getattr(obj, atr), 'error', error_dict[label][atr])
            if label in fixed_dict.keys():
                for atr in fixed_dict[label].keys():
                    setattr(getattr(obj, atr), 'fixed', fixed_dict[label][atr])
            atoms.append(obj)
        return atoms

    def from_cif_block(self, block: gemmi.cif.Block) -> B:
        atoms = self._site_runner(block)
        return self._CIF_CLASS('from_cif', *atoms)

    def add_to_cif_block(self, obj: B, block: gemmi.cif.Block) -> NoReturn:
        additional_keys = []
        additional_objs = []

        ADP_WRITER = AtomicDisplacement()
        if len(obj) > 0:
            if getattr(obj[0], 'adp', False):
                # additional keys for isotropic adp
                additional_keys, additional_objs = ADP_WRITER.iso_adp_block(obj, block)
        # Add main atom loop
        self._add_to_cif_block(obj, block, additional_keys, additional_objs)
        # another loop for anisotropic atomic displacements
        ADP_WRITER.add_to_cif_block(obj, block)
        # another loop for anisotropic magnetic susceptibilities
        MSP_WRITER = MagneticSusceptibility()
        MSP_WRITER.add_to_cif_block(obj, block)

    def _add_to_cif_block(self, obj: B, block: gemmi.cif.Block, additional_keys, additional_objs) -> NoReturn:
        # First add the main loop
        items = list(self._CIF_CONVERSIONS)
        names = [item[1] for item in items]
        lines = []
        # this loop adds atom name, symbol and coordinates as a line under the main "loop_"
        for idx1, atom in enumerate(obj):
            line = []
            for idx2, item in enumerate(items):
                if item[0] == 'specie':
                    s = str(atom.specie)
                else:
                    s = self.variable_to_string(atom.__getattribute__(item[0]))
                line.append(s)
            lines.append(line)

        # this loop adds ADP values to the end of the atom description lines
        for keys, objs, line in zip(additional_keys, additional_objs, lines):
            for idx, _obj in enumerate(objs):
                for key in keys[idx]:
                    k = key[0]
                    # Need to assure Uiso is present on ADP, else just assume 0.0
                    s = self.variable_to_string(_obj.__getattribute__(k)) if hasattr(_obj, k) else '0.0'
                    # for the cif representation, we only use [UB]iso
                    if 'ani' in s:
                        s = s.replace('ani', 'iso')
                    line.append(s)
                    if key[1] not in names:
                        names.append(key[1])

        # init_loop creates the main "loop_" section with '_atom_site_<names>' as the tags
        loop = block.init_loop(self._CIF_SECTION_NAME, names)
        for line in lines:
            loop.add_row(line)

    def from_cif_string(self, cif_string: str) -> List[B]:
        if 'data_' not in cif_string:
            cif_string = 'data_temp\n' + cif_string
        cif_blocks = gemmi.cif.read_string(cif_string)
        objs = []
        for block in cif_blocks:
            objs.append(self.from_cif_block(block))
        return objs
