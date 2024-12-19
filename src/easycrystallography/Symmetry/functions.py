ACCESIBLE_IT_NUMBER_TRICLINIC_SYSTEM = tuple(range(1, 3))
ACCESIBLE_IT_NUMBER_MONOCLINIC_SYSTEM = tuple(range(3, 16))
ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM = tuple(range(16, 75))
ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM = tuple(range(7, 143))
ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM = tuple(range(143, 168))
ACCESIBLE_IT_NUMBER_HEXAGONAL_SYSTEM = tuple(range(168, 195))
ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM = tuple(range(195, 231))

ACCESIBLE_IT_NUMBER_MONOCLINIC_SYSTEM_TRIPLE_CHOICE = (5, 7, 8, 9, 12, 13, 14, 15)
ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM_DOUBLE_CHOICE = (48, 50, 59, 68, 70)
ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM_DOUBLE_CHOICE = (85, 86, 88, 125, 126, 129, 130, 133, 134, 137, 138, 141, 142)
ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM_DOUBLE_AXES = (146, 148, 155, 160, 161, 166, 167)
ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM_DOUBLE_CHOICE = (201, 203, 222, 224, 227, 228)

ACCESIBLE_IT_NUMBER = (
    ACCESIBLE_IT_NUMBER_TRICLINIC_SYSTEM
    + ACCESIBLE_IT_NUMBER_MONOCLINIC_SYSTEM
    + ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM
    + ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM
    + ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM
    + ACCESIBLE_IT_NUMBER_HEXAGONAL_SYSTEM
    + ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM
)


def get_default_it_coordinate_system_code_by_it_number(it_number: int) -> str:
    crystal_system = get_crystal_system_by_it_number(it_number)
    if crystal_system == 'triclinic':
        it_coordinate_system_code = None
    elif crystal_system == 'monoclinic':
        it_coordinate_system_code = 'b1'
    elif crystal_system == 'orthorhombic':
        if it_number in ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_code = '2abc'
        else:
            it_coordinate_system_code = 'abc'
    elif crystal_system == 'tetragonal':
        if it_number in ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_code = '2'
        else:
            it_coordinate_system_code = '1'
    elif crystal_system == 'trigonal':
        if it_number in ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM_DOUBLE_AXES:
            it_coordinate_system_code = 'h'
        else:
            it_coordinate_system_code = 'r'
    elif crystal_system == 'hexagonal':
        it_coordinate_system_code = 'h'
    elif crystal_system == 'cubic':
        if it_number in ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_code = '2'
        else:
            it_coordinate_system_code = '1'
    else:
        it_coordinate_system_code = None
    return it_coordinate_system_code


def get_crystal_system_by_it_number(it_number: int) -> str:
    if it_number is None:
        return None
    if (it_number >= 1) & (it_number <= 2):
        res = 'triclinic'
    elif (it_number >= 3) & (it_number <= 15):
        res = 'monoclinic'
    elif (it_number >= 16) & (it_number <= 74):
        res = 'orthorhombic'
    elif (it_number >= 75) & (it_number <= 142):
        res = 'tetragonal'
    elif (it_number >= 143) & (it_number <= 167):
        res = 'trigonal'
    elif (it_number >= 168) & (it_number <= 194):
        res = 'hexagonal'
    elif (it_number >= 195) & (it_number <= 230):
        res = 'cubic'
    else:
        res = None
    return res
