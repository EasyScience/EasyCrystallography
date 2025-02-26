# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

import re
from copy import deepcopy
from typing import Any
from typing import Dict
from typing import Union

import periodictable as pt
from easyscience.Objects.variable import DescriptorStr as Descriptor
from easyscience.Utils.classTools import addProp
from easyscience.Utils.classTools import removeProp

_SPECIE_DETAILS = {
    'type_symbol': {
        'description': 'A code to identify the atom species occupying this site.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_atom_site/',
    },
}

_REDIRECT = deepcopy(Descriptor._REDIRECT)
_REDIRECT['specie'] = lambda obj: obj._raw_data['str']


class Specie(Descriptor):
    _REDIRECT = _REDIRECT

    def __init__(self, specie: str = 'H', **kwargs):
        if 'value' in kwargs.keys():
            specie = kwargs.pop('value')

        self._raw_data: Dict[str, Union[str, int]] = {}
        self._props: Dict[str, Any] = {}

        self._reset_data()
        super(Specie, self).__init__('specie', self.__gen_data(specie), **_SPECIE_DETAILS['type_symbol'])

        # Monkey patch the unit and the value to take into account the new type situation
        self.__previous_set = self.__class__.value.fset

        addProp(
            self,
            'value',
            fget=self.__class__.value.fget,
            fset=lambda obj, val: self.__previous_set(obj, obj.__gen_data(val)),
            fdel=self.__class__.value.fdel,
        )

    def _reset_data(self):
        self._raw_data = dict.fromkeys(['str', 'observed', 'element', 'isotope', 'oxi_state', 'spin'])
        self._props = {}

    def __gen_data(self, value_str: str):
        s = re.search(r'([0-9.]*)([A-Z][a-z]*)([0-9.]*)([+\-]*)', value_str)
        # group(1) = Isotope
        # group(2) = Element
        # group(3) = Oxi state
        # group(4) = Oxi +/2

        if s is None:
            raise ValueError('Invalid specie string')
        element_str = s.group(2)
        element = getattr(pt, element_str, None)
        if element is None:
            raise ValueError(f'Element ({s.group(2)}) not found in periodictable')

        if self._props is not None:
            for k in self._props.keys():
                removeProp(self, k)

        self._reset_data()
        self._raw_data['str'] = value_str
        self._raw_data['element'] = element
        self._raw_data['observed'] = element

        props = {
            s: getattr(self._raw_data['observed'], s)
            for s in self._raw_data['observed'].__dir__()
            if not s.startswith('_') and hasattr(self._raw_data['observed'], s)
        }
        isotope_str = s.group(1)
        self._raw_data['isotope'] = None
        if isotope_str:
            self._raw_data['isotope'] = int(isotope_str)
            self._raw_data['observed'] = pt.core.Isotope(self._raw_data['observed'], self._raw_data['isotope'])
            props.update(
                {
                    s: getattr(self._raw_data['observed'], s)
                    for s in self._raw_data['observed'].__dir__()
                    if not s.startswith('_') and hasattr(self._raw_data['observed'], s)
                }
            )
        oxi_state_str = s.group(3)
        oxi_pm_str = s.group(4)
        self._raw_data['oxi_state'] = None
        if oxi_state_str:
            self._raw_data['oxi_state'] = int(oxi_pm_str + oxi_state_str)
            self._raw_data['observed'] = pt.core.Ion(self._raw_data['observed'], self._raw_data['oxi_state'])
            props.update(
                {
                    s: getattr(self._raw_data['observed'], s)
                    for s in self._raw_data['observed'].__dir__()
                    if not s.startswith('_') and hasattr(self._raw_data['observed'], s)
                }
            )
        props['common_name'] = props.pop('name')
        self._props = props
        for k, v in props.items():
            addProp(self, k, fget=self.__getter_periodic(k))
        rep = f'{self._raw_data["element"]}'
        if self.is_ion:
            rep += f'{abs(self._raw_data["oxi_state"])}'
            if self._raw_data['oxi_state'] > 0:
                rep += '+'
            else:
                rep += '-'
        return rep

    def __repr__(self) -> str:
        rep = f"<{self.__class__.__name__} '{self.name}': "
        if self.is_isotope:
            rep += f'{self._raw_data["isotope"]}'
        rep += f'{self._raw_data["element"]}'
        if self.is_ion:
            if self._raw_data['oxi_state'] > 0:
                rep += f'{self._raw_data["oxi_state"]}+'
            else:
                rep += f'{self._raw_data["oxi_state"]}-'
        rep += '>'
        return rep

    def __str__(self) -> str:
        rep = ''
        if self.is_isotope:
            rep += f'{self._raw_data["isotope"]}'
        rep += f'{self._raw_data["element"]}'
        if self.is_ion:
            if self._raw_data['oxi_state'] > 0:
                rep += f'{self._raw_data["oxi_state"]}+'
            else:
                rep += f'{str(self._raw_data["oxi_state"])[1:]}-'
        return rep

    @staticmethod
    def __getter_periodic(key: str):
        def getter(obj):
            return obj._props.get(key)

        return getter

    @property
    def is_ion(self):
        return self._raw_data['oxi_state'] is not None

    @property
    def is_isotope(self):
        return self._raw_data['isotope'] is not None
