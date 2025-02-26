# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

from __future__ import annotations

from typing import TYPE_CHECKING
from typing import ClassVar
from typing import List
from typing import Optional
from typing import Type
from typing import Union

import numpy as np
from easyscience.Objects.ObjectClasses import BaseObj
from easyscience.Objects.variable import DescriptorStr as Descriptor
from easyscience.Objects.variable import Parameter
from easyscience.Utils.classTools import addProp
from easyscience.Utils.classTools import removeProp

if TYPE_CHECKING:
    from easyscience.Utils.typing import iF

_ANIO_DETAILS = {
    'msp_type': {
        'description': 'A standard code used to describe the type of atomic displacement parameters used for the site.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_atom_site/',
        'value': 'Uani',
    },
    'Cani': {
        'description': 'Isotropic magnetic susceptibility parameter.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_atom_site/',
        'value': 0.0,
        'unit': 'T^-1',
        'fixed': True,
    },
    'Ciso': {
        'description': 'Isotropic magnetic susceptibility parameter, or equivalent isotropic magnetic susceptibility '
        'parameter, C(equiv), in inverted teslas, calculated from anisotropic susceptibility '
        'components.',
        'value': 0.0,
        'max': np.inf,
        'unit': 'T^-1',
        'fixed': True,
    },
}


class MSPBase(BaseObj):
    def __init__(self, *args, **kwargs):
        super(MSPBase, self).__init__(*args, **kwargs)

    @property
    def matrix(self) -> np.ndarray:
        matrix = np.zeros([3, 3])
        pars = self.get_parameters()
        if len(pars) == 1:
            np.fill_diagonal(matrix, pars[0].value)
        elif len(pars) == 6:
            matrix[0, 0] = pars[0].value
            matrix[0, 1] = pars[1].value
            matrix[0, 2] = pars[2].value
            matrix[1, 1] = pars[3].value
            matrix[1, 2] = pars[4].value
            matrix[2, 2] = pars[5].value
        return matrix


class Cani(MSPBase):
    chi_11: ClassVar[Parameter]
    chi_12: ClassVar[Parameter]
    chi_13: ClassVar[Parameter]
    chi_22: ClassVar[Parameter]
    chi_23: ClassVar[Parameter]
    chi_33: ClassVar[Parameter]

    def __init__(
        self,
        chi_11: Optional[Union[Parameter, float]] = None,
        chi_12: Optional[Union[Parameter, float]] = None,
        chi_13: Optional[Union[Parameter, float]] = None,
        chi_22: Optional[Union[Parameter, float]] = None,
        chi_23: Optional[Union[Parameter, float]] = None,
        chi_33: Optional[Union[Parameter, float]] = None,
        msp_values: Optional[Type[MSPBase]] = None,
        interface=None,
    ):
        super(Cani, self).__init__(
            'Cani',
            chi_11=Parameter('chi_11', **_ANIO_DETAILS['Cani']),
            chi_12=Parameter('chi_12', **_ANIO_DETAILS['Cani']),
            chi_13=Parameter('chi_13', **_ANIO_DETAILS['Cani']),
            chi_22=Parameter('chi_22', **_ANIO_DETAILS['Cani']),
            chi_23=Parameter('chi_23', **_ANIO_DETAILS['Cani']),
            chi_33=Parameter('chi_33', **_ANIO_DETAILS['Cani']),
        )
        if chi_11 is not None:
            self.chi_11 = chi_11
        if chi_12 is not None:
            self.chi_12 = chi_12
        if chi_13 is not None:
            self.chi_13 = chi_13
        if chi_22 is not None:
            self.chi_22 = chi_22
        if chi_23 is not None:
            self.chi_23 = chi_23
        if chi_33 is not None:
            self.chi_33 = chi_33
        if msp_values is not None:
            self.chi_11 = msp_values.chi_11
            self.chi_12 = msp_values.chi_12
            self.chi_13 = msp_values.chi_13
            self.chi_22 = msp_values.chi_22
            self.chi_23 = msp_values.chi_23
            self.chi_33 = msp_values.chi_33
        self.interface = interface


class Ciso(MSPBase):
    chi: ClassVar[Parameter]

    def __init__(
        self,
        chi: Optional[Union[Parameter, float]] = None,
        msp_values: Optional[Type[MSPBase]] = None,
        interface: Optional[iF] = None,
    ):
        super(Ciso, self).__init__('Ciso', chi=Parameter('chi', **_ANIO_DETAILS['Ciso']))
        if chi is not None:
            self.chi = chi
        if msp_values is not None:
            self.chi = msp_values.chi
        self.interface = interface


_AVAILABLE_ISO_TYPES = {'Cani': Cani, 'Ciso': Ciso}


class MagneticSusceptibility(BaseObj):
    msp_type: ClassVar[Descriptor]
    msp_class: ClassVar[Type[MSPBase]]

    def __init__(self, msp_type: Union[Descriptor, str], interface: Optional = None, **kwargs):
        if isinstance(msp_type, str):
            msp_type = Descriptor('msp_type', msp_type)
        msp_class_name = msp_type.value
        if msp_class_name in _AVAILABLE_ISO_TYPES.keys():
            msp_class = _AVAILABLE_ISO_TYPES[msp_class_name]
            if 'msp_class' in kwargs:
                # enable passing chi values directly to constructor
                kwargs['msp_values'] = kwargs['msp_class']
                _ = kwargs.pop('msp_class')
            msp = msp_class(**kwargs, interface=interface)
        else:
            raise AttributeError(f'{msp_class_name} is not a valid magnetic susceptibility type')
        super(MagneticSusceptibility, self).__init__('msp', msp_type=msp_type, msp_class=msp)
        for par in msp.get_parameters():
            addProp(
                self,
                par.name,
                fget=self.__a_getter(par.name),
                fset=self.__a_setter(par.name),
            )
        self.interface = interface

    def switch_type(self, msp_string: str, **kwargs):
        if msp_string in _AVAILABLE_ISO_TYPES.keys():
            msp_class = _AVAILABLE_ISO_TYPES[msp_string]
            if kwargs:
                msp_class: MSPBase = msp_class(interface=self.interface, **kwargs)
            else:
                msp_class: MSPBase = msp_class(interface=self.interface)
        else:
            raise AttributeError

        for par in self.msp_class.get_parameters():
            removeProp(self, par.name)
        self.msp_class = msp_class
        self.msp_type = msp_string
        for par in msp_class.get_parameters():
            addProp(
                self,
                par.name,
                fget=self.__a_getter(par.name),
                fset=self.__a_setter(par.name),
            )

    @property
    def available_types(self) -> List[str]:
        return [name for name in _AVAILABLE_ISO_TYPES.keys()]

    @staticmethod
    def __a_getter(key: str):
        def getter(obj):
            return obj.msp_class._kwargs[key]

        return getter

    @staticmethod
    def __a_setter(key):
        def setter(obj, value):
            obj.msp_class._kwargs[key].value = value

        return setter
