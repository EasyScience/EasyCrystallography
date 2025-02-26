# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

from __future__ import annotations

import textwrap
from abc import abstractmethod
from math import floor
from math import log10
from numbers import Number
from typing import TYPE_CHECKING
from typing import List
from typing import NoReturn
from typing import Optional

import gemmi
import numpy as np

if TYPE_CHECKING:
    from easyscience.Utils.typing import B
    from easyscience.Utils.typing import V

_MAX_LEN = 140
_MAX_LABEL_LEN = 130


class CIF_Template:
    """
    This class is an abstract class for a CIF template.
    """

    def __init__(self, decimal_places: int = 8):
        self.decimal_places = decimal_places

    @property
    def CLASS_READER(self):
        return self.from_cif_block

    @property
    def CLASS_WRITER(self):
        return self.add_to_cif_block

    @property
    def STRING_READER(self):
        return self.from_cif_string

    @property
    def STRING_WRITER(self):
        return self.to_cif_str

    @abstractmethod
    def from_cif_block(self, block: gemmi.cif.Block) -> B:
        pass

    def to_cif_str(self, obj: B, block: Optional[gemmi.cif.Block] = None) -> str:
        """
        This method returns the CIF string.
        :return: The CIF string.
        """
        if block is None:
            block = gemmi.cif.Block('temp')
        start_str = block.as_string()
        self.add_to_cif_block(obj, block)
        final_str = block.as_string()
        return final_str[len(start_str) :]

    @abstractmethod
    def add_to_cif_block(self, obj: B, block: gemmi.cif.Block) -> NoReturn:
        pass

    @abstractmethod
    def from_cif_string(self, cif_string: str) -> List[B]:
        pass

    def variable_to_string(self, obj: V) -> str:
        value = obj.value
        error = getattr(obj, 'error', None)
        if error is not None and np.isclose(error, 0.0):
            error = None
        decimal_places = self.decimal_places
        if isinstance(obj.value, Number):
            if abs(value - int(value)) > 0.0:
                decimal_places = len(str(value).split('.')[1])
            initial_str = '{:.' + str(decimal_places) + 'f}'
            s = initial_str.format(round(value, decimal_places))
            if error is not None:
                xe_exp = int(floor(log10(error)))

                # uncertainty
                un_exp = xe_exp - decimal_places + 1
                un_int = round(error * 10 ** (-un_exp))

                # nominal value
                no_exp = un_exp
                no_int = round(value * 10 ** (-no_exp))

                # format - nom(unc)
                fmt = '%%.%df' % max(0, -no_exp)
                s = (fmt + '(%.0f)') % (
                    no_int * 10**no_exp,
                    un_int * 10 ** max(0, un_exp),
                )
        elif isinstance(value, str):
            s = '{:s}'.format(value)
        else:
            s = '{:s}'.format(str(value))
        if getattr(obj, 'fixed', None) is not None and not obj.fixed and error is None:
            s += '()'
        return CIF_Template._format_field(s)

    @staticmethod
    def string_to_variable(in_string: str):
        in_string = in_string.strip()
        if "'" in in_string:
            in_string = in_string.replace("'", '')
        if '"' in in_string:
            in_string = in_string.replace('"', '')
        fixed = None
        error = None
        tokens = in_string.split('(')
        try:
            value = float(tokens[0])
        except ValueError:
            value = tokens[0]
            return value, error, fixed
        if len(tokens) > 1:
            fixed = False
            if tokens[1][0] != ')':
                decimal_places = len(tokens[0].split('.')[1])
                error = int(tokens[1][:-1]) * (10**-decimal_places)
                # error = (10 ** -(len(f'{tokens[0]}'.split('.')[1]) - 1)) * int(tokens[1][:-1])
        return value, error, fixed

    @staticmethod
    def _format_field(v):
        v = v.__str__().strip()
        if len(v) > _MAX_LEN:
            return ';\n' + textwrap.fill(v, _MAX_LEN) + '\n;'
        # add quotes if necessary
        if v == '':
            return '""'
        if (' ' in v or v[0] == '_') and not (v[0] == "'" and v[-1] == "'") and not (v[0] == '"' and v[-1] == '"'):
            if "'" in v:
                q = '"'
            else:
                q = "'"
            v = q + v + q
        return v
