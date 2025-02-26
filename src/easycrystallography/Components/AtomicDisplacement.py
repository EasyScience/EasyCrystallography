# SPDX-FileCopyrightText: 2024 EasyCrystallography contributors
# SPDX-License-Identifier: BSD-3-Clause
# © 2022-2024 Contributors to the EasyCrystallography project <https://github.com/EasyScience/EasyCrystallography>

from __future__ import annotations

from typing import TYPE_CHECKING
from typing import ClassVar
from typing import List
from typing import Optional
from typing import TypeVar
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
    'adp_type': {
        'description': 'A standard code used to describe the type of atomic displacement parameters used for the site.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_atom_site/',
        'value': 'Uani',
    },
    'Uani': {
        'description': 'Isotropic atomic displacement parameter, or equivalent isotropic atomic  displacement '
        'parameter, U(equiv), in angstroms squared, calculated from anisotropic atomic displacement  '
        'parameters.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_atom_site/',
        'value': 0.0,
        'unit': 'angstrom^2',
        'fixed': True,
    },
    'Uiso': {
        'description': 'The standard anisotropic atomic displacement components in angstroms squared which appear in '
        'the structure-factor term.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_atom_site/',
        'value': 0.0,
        'min': 0,
        'max': np.inf,
        'unit': 'angstrom^2',
        'fixed': True,
    },
    'Bani': {
        'description': 'The standard anisotropic atomic displacement components in angstroms squared which appear in '
        'the structure-factor term.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_atom_site/',
        'value': 0.0,
        'unit': 'angstrom^2',
        'fixed': True,
    },
    'Biso': {
        'description': 'Isotropic atomic displacement parameter, or equivalent isotropic atomic displacement '
        'parameter, B(equiv), in angstroms squared, calculated from anisotropic displacement '
        'components.',
        'url': 'https://docs.easydiffraction.org/lib/dictionaries/_atom_site/',
        'value': 0.0,
        'min': 0,
        'max': np.inf,
        'unit': 'angstrom^2',
        'fixed': True,
    },
}


class AdpBase(BaseObj):
    def __init__(self, *args, **kwargs):
        super(AdpBase, self).__init__(*args, **kwargs)

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

    def __repr__(self):
        s = f'{self.name} - ('
        for par in self.get_parameters():
            s += f'{par.name}: {par.value}, '
        s = s[:-2] + ')'
        return s


class Anisotropic(AdpBase):
    U_11: ClassVar[Parameter]
    U_12: ClassVar[Parameter]
    U_13: ClassVar[Parameter]
    U_22: ClassVar[Parameter]
    U_23: ClassVar[Parameter]
    U_33: ClassVar[Parameter]

    def __init__(
        self,
        Uiso_ani: Optional[float] = None,
        U_11: Optional[Union[Parameter, float]] = None,
        U_12: Optional[Union[Parameter, float]] = None,
        U_13: Optional[Union[Parameter, float]] = None,
        U_22: Optional[Union[Parameter, float]] = None,
        U_23: Optional[Union[Parameter, float]] = None,
        U_33: Optional[Union[Parameter, float]] = None,
        interface: Optional[iF] = None,
    ):
        super(Anisotropic, self).__init__(
            'anisoU',
            U_11=Parameter('U_11', **_ANIO_DETAILS['Uani']),
            U_12=Parameter('U_12', **_ANIO_DETAILS['Uani']),
            U_13=Parameter('U_13', **_ANIO_DETAILS['Uani']),
            U_22=Parameter('U_22', **_ANIO_DETAILS['Uani']),
            U_23=Parameter('U_23', **_ANIO_DETAILS['Uani']),
            U_33=Parameter('U_33', **_ANIO_DETAILS['Uani']),
        )
        if U_11 is not None:
            self.U_11 = U_11
        if U_12 is not None:
            self.U_12 = U_12
        if U_13 is not None:
            self.U_13 = U_13
        if U_22 is not None:
            self.U_22 = U_22
        if U_23 is not None:
            self.U_23 = U_23
        if U_33 is not None:
            self.U_33 = U_33
        if Uiso_ani is not None:
            self.Uiso_ani = Uiso_ani
        else:
            # for cubic, tetragonal, and orthorhombic systems
            self.Uiso_ani = (self.U_11.value + self.U_22.value + self.U_33.value) / 3.0
        self.interface = interface


class Isotropic(AdpBase):
    Uiso: ClassVar[Parameter]

    def __init__(
        self,
        Uiso: Optional[Union[Parameter, float]] = None,
        interface: Optional[iF] = None,
    ):
        super(Isotropic, self).__init__('Uiso', Uiso=Parameter('Uiso', **_ANIO_DETAILS['Uiso']))
        if Uiso is not None:
            self.Uiso = Uiso
        self.interface = interface


class AnisotropicBij(AdpBase):
    B_11: ClassVar[Parameter]
    B_12: ClassVar[Parameter]
    B_13: ClassVar[Parameter]
    B_22: ClassVar[Parameter]
    B_23: ClassVar[Parameter]
    B_33: ClassVar[Parameter]

    def __init__(
        self,
        B_11: Optional[Union[Parameter, float]] = None,
        B_12: Optional[Union[Parameter, float]] = None,
        B_13: Optional[Union[Parameter, float]] = None,
        B_22: Optional[Union[Parameter, float]] = None,
        B_23: Optional[Union[Parameter, float]] = None,
        B_33: Optional[Union[Parameter, float]] = None,
        interface: Optional[iF] = None,
    ):
        super(AnisotropicBij, self).__init__(
            'anisoB',
            **{name: Parameter(name, **_ANIO_DETAILS['Bani']) for name in ['B_11', 'B_12', 'B_13', 'B_22', 'B_23', 'B_33']},
        )
        if B_11 is not None:
            self.B_11 = B_11
        if B_12 is not None:
            self.B_12 = B_12
        if B_13 is not None:
            self.B_13 = B_13
        if B_22 is not None:
            self.B_22 = B_22
        if B_23 is not None:
            self.B_23 = B_23
        if B_33 is not None:
            self.B_33 = B_33
        self.interface = interface


class IsotropicB(AdpBase):
    Biso: ClassVar[Parameter]

    def __init__(
        self,
        Biso: Optional[Union[Parameter, float]] = None,
        interface: Optional[iF] = None,
    ):
        super(IsotropicB, self).__init__('Biso', Biso=Parameter('Biso', **_ANIO_DETAILS['Biso']))
        if Biso is not None:
            self.Biso = Biso
        self.interface = interface


_AVAILABLE_ISO_TYPES = {
    'Uani': Anisotropic,
    'Uiso': Isotropic,
    # 'Uovl': 'Overall',
    # 'Umpe': 'MultipoleExpansion',
    'Bani': AnisotropicBij,
    'Biso': IsotropicB,
    # 'Bovl': 'OverallB'
}

if TYPE_CHECKING:
    AB = TypeVar('AB', bound=AdpBase)


class AtomicDisplacement(BaseObj):
    adp_type: ClassVar[Descriptor]
    adp_class: ClassVar[AB]

    def __init__(
        self,
        adp_type: Optional[Union[Descriptor, str]] = None,
        interface: Optional[iF] = None,
        **kwargs,
    ):
        if adp_type is None:
            adp_type = 'Uiso'
        if isinstance(adp_type, str):
            adp_type = Descriptor('adp_type', adp_type)
        adp_class_name = adp_type.value
        if adp_class_name in _AVAILABLE_ISO_TYPES.keys():
            adp_class = _AVAILABLE_ISO_TYPES[adp_class_name]
            # enable passing ADP parameters to constructor
            if 'adp_class' in kwargs.keys():
                m = getattr(kwargs['adp_class'], adp_class_name)
                kwargs[adp_class_name] = m
                _ = kwargs.pop('adp_class')
            adp = adp_class(**kwargs, interface=interface)
        else:
            raise AttributeError(f'{adp_class_name} is not a valid adp type')
        super(AtomicDisplacement, self).__init__('adp', adp_type=adp_type, adp_class=adp)
        for par in adp.get_parameters():
            addProp(
                self,
                par.name,
                fget=self.__a_getter(par.name),
                fset=self.__a_setter(par.name),
            )
        self.interface = interface

    def switch_type(self, adp_string: str, **kwargs):
        # if adp_string in _AVAILABLE_ISO_TYPES.keys():
        # adp_class = _AVAILABLE_ISO_TYPES[adp_string]
        # if kwargs:
        #     adp = adp_class(**kwargs, interface=self.interface)
        # else:
        #     adp = adp_class(interface=self.interface)
        # else:
        #     raise AttributeError(f"{adp_string} is not a valid adp type")
        for par in self.adp_class.get_parameters():
            removeProp(self, par.name)
        self.__init__(adp_type=adp_string, **kwargs)

    @property
    def available_types(self) -> List[str]:
        return [name for name in _AVAILABLE_ISO_TYPES.keys()]

    @staticmethod
    def __a_getter(key: str):
        def getter(obj):
            return obj.adp_class._kwargs[key]

        return getter

    @staticmethod
    def __a_setter(key):
        def setter(obj, value):
            obj.adp_class._kwargs[key].value = value

        return setter
