# ruff: noqa
__author__ = 'Shyue Ping Ong, Shyam Dwaraknath, Matthew Horton'
__copyright__ = 'Copyright 2011, The Materials Project'
__version__ = '1.0'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'shyuep@gmail.com'
__status__ = 'Production'
__date__ = 'Sep 23, 2011'

from typing import List
from typing import Tuple
from typing import Union

"""
This module provides classes that operate on points or vectors in 3D space.
"""

import re
import string
import warnings
from math import cos
from math import pi
from math import sin
from math import sqrt

import numpy as np
from easyscience.Objects.core import ComponentSerializer
from easyscience.Utils.string import transformation_to_string


class SymmOp(ComponentSerializer):
    """
    A symmetry operation in cartesian space. Consists of a rotation plus a
    translation. Implementation is as an affine transformation matrix of rank 4
    for efficiency. Read: http://en.wikipedia.org/wiki/Affine_transformation.
    .. attribute:: affine_matrix
        A 4x4 numpy.array representing the symmetry operation.
    """

    _REDIRECT = {'affine_transformation_matrix': lambda obj: getattr(obj, 'affine_matrix').aslist()}

    def __init__(self, affine_transformation_matrix, tol=0.01):
        """
        Initializes the SymmOp from a 4x4 affine transformation matrix.
        In general, this constructor should not be used unless you are
        transferring rotations.  Use the static constructors instead to
        generate a SymmOp from proper rotations and translation.
        Args:
            affine_transformation_matrix (4x4 array): Representing an
                affine transformation.
            tol (float): Tolerance for determining if matrices are equal.
        """
        affine_transformation_matrix = np.array(affine_transformation_matrix)
        if affine_transformation_matrix.shape != (4, 4):
            raise ValueError('Affine Matrix must be a 4x4 numpy array!')
        self.affine_matrix = affine_transformation_matrix
        self.tol = tol

    @staticmethod
    def from_rotation_and_translation(
        rotation_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        translation_vec=(0, 0, 0),
        tol=0.1,
    ):
        """
        Creates a symmetry operation from a rotation matrix and a translation
        vector.
        Args:
            rotation_matrix (3x3 array): Rotation matrix.
            translation_vec (3x1 array): Translation vector.
            tol (float): Tolerance to determine if rotation matrix is valid.
        Returns:
            SymmOp object
        """
        rotation_matrix = np.array(rotation_matrix)
        translation_vec = np.array(translation_vec)
        if rotation_matrix.shape != (3, 3):
            raise ValueError('Rotation Matrix must be a 3x3 numpy array.')
        if translation_vec.shape != (3,):
            raise ValueError('Translation vector must be a rank 1 numpy array with 3 elements.')
        affine_matrix = np.eye(4)
        affine_matrix[0:3][:, 0:3] = rotation_matrix
        affine_matrix[0:3][:, 3] = translation_vec
        return SymmOp(affine_matrix, tol)

    def __eq__(self, other):
        return np.allclose(self.affine_matrix, other.affine_matrix, atol=self.tol)

    def __hash__(self):
        return 7

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        output = [
            'Rot:',
            str(self.affine_matrix[0:3][:, 0:3]),
            'tau',
            str(self.affine_matrix[0:3][:, 3]),
        ]
        return '\n'.join(output)

    def operate(self, point):
        """
        Apply the operation on a point.
        Args:
            point: Cartesian coordinate.
        Returns:
            Coordinates of point after operation.
        """
        affine_point = np.array([point[0], point[1], point[2], 1])
        return np.dot(self.affine_matrix, affine_point)[0:3]

    def operate_multi(self, points):
        """
        Apply the operation on a list of points.
        Args:
            points: List of Cartesian coordinates
        Returns:
            Numpy array of coordinates after operation
        """
        points = np.array(points)
        affine_points = np.concatenate([points, np.ones(points.shape[:-1] + (1,))], axis=-1)
        return np.inner(affine_points, self.affine_matrix)[..., :-1]

    def apply_rotation_only(self, vector):
        """
        Vectors should only be operated by the rotation matrix and not the
        translation vector.
        Args:
            vector (3x1 array): A vector.
        """
        return np.dot(self.rotation_matrix, vector)

    def transform_tensor(self, tensor):
        """
        Applies rotation portion to a tensor. Note that tensor has to be in
        full form, not the Voigt form.
        Args:
            tensor (numpy array): a rank n tensor
        Returns:
            Transformed tensor.
        """
        dim = tensor.shape
        rank = len(dim)
        if not all([i == 3 for i in dim]):
            raise ValueError('All dimensions must be equal to 3')

        # Build einstein sum string
        lc = string.ascii_lowercase
        indices = lc[:rank], lc[rank : 2 * rank]
        einsum_string = ','.join([a + i for a, i in zip(*indices)])
        einsum_string += ',{}->{}'.format(*indices[::-1])
        einsum_args = [self.rotation_matrix] * rank + [tensor]

        return np.einsum(einsum_string, *einsum_args)

    def are_symmetrically_related(self, point_a, point_b, tol=0.001):
        """
        Checks if two points are symmetrically related.
        Args:
            point_a (3x1 array): First point.
            point_b (3x1 array): Second point.
            tol (float): Absolute tolerance for checking distance.
        Returns:
            True if self.operate(point_a) == point_b or vice versa.
        """
        if np.allclose(self.operate(point_a), point_b, atol=tol):
            return True
        if np.allclose(self.operate(point_b), point_a, atol=tol):
            return True
        return False

    @property
    def rotation_matrix(self) -> np.ndarray:
        """
        A 3x3 numpy.array representing the rotation matrix.
        """
        return self.affine_matrix[0:3][:, 0:3]

    @property
    def translation_vector(self):
        """
        A rank 1 numpy.array of dim 3 representing the translation vector.
        """
        return self.affine_matrix[0:3][:, 3]

    def __mul__(self, other):
        """
        Returns a new SymmOp which is equivalent to apply the "other" SymmOp
        followed by this one.
        """
        new_matrix = np.dot(self.affine_matrix, other.affine_matrix)
        return SymmOp(new_matrix)

    @property
    def inverse(self):
        """
        Returns inverse of transformation.
        """
        invr = np.linalg.inv(self.affine_matrix)
        return SymmOp(invr)

    @staticmethod
    def from_axis_angle_and_translation(
        axis: Union[Tuple[float, float, float], List[float], np.ndarray],
        angle: float,
        angle_in_radians: bool = False,
        translation_vec: Union[Tuple[float, float, float], List[float], np.ndarray] = (
            0,
            0,
            0,
        ),
    ):
        """
        Generates a SymmOp for a rotation about a given axis plus translation.
        Args:
            axis: The axis of rotation in cartesian space. For example,
                [1, 0, 0] indicates rotation about x-axis.
            angle (float): Angle of rotation.
            angle_in_radians (bool): Set to True if angles are given in
                radians. Or else, units of degrees are assumed.
            translation_vec: A translation vector. Defaults to zero.
        Returns:
            SymmOp for a rotation about given axis and translation.
        """
        if isinstance(axis, (tuple, list)):
            axis = np.array(axis)

        if isinstance(translation_vec, (tuple, list)):
            vec = np.array(translation_vec)
        else:
            vec = translation_vec

        a = angle if angle_in_radians else angle * pi / 180
        cosa = cos(a)
        sina = sin(a)
        u = axis / np.linalg.norm(axis)
        r = np.zeros((3, 3))
        r[0, 0] = cosa + u[0] ** 2 * (1 - cosa)
        r[0, 1] = u[0] * u[1] * (1 - cosa) - u[2] * sina
        r[0, 2] = u[0] * u[2] * (1 - cosa) + u[1] * sina
        r[1, 0] = u[0] * u[1] * (1 - cosa) + u[2] * sina
        r[1, 1] = cosa + u[1] ** 2 * (1 - cosa)
        r[1, 2] = u[1] * u[2] * (1 - cosa) - u[0] * sina
        r[2, 0] = u[0] * u[2] * (1 - cosa) - u[1] * sina
        r[2, 1] = u[1] * u[2] * (1 - cosa) + u[0] * sina
        r[2, 2] = cosa + u[2] ** 2 * (1 - cosa)

        return SymmOp.from_rotation_and_translation(r, vec)

    @staticmethod
    def from_origin_axis_angle(
        origin: Union[Tuple[float, float, float], List[float], np.ndarray],
        axis: Union[Tuple[float, float, float], List[float], np.ndarray],
        angle: float,
        angle_in_radians: bool = False,
    ):
        """
        Generates a SymmOp for a rotation about a given axis through an
        origin.
        Args:
            origin (3x1 array): The origin which the axis passes through.
            axis (3x1 array): The axis of rotation in cartesian space. For
                example, [1, 0, 0]indicates rotation about x-axis.
            angle (float): Angle of rotation.
            angle_in_radians (bool): Set to True if angles are given in
                radians. Or else, units of degrees are assumed.
        Returns:
            SymmOp.
        """
        theta = angle * pi / 180 if not angle_in_radians else angle
        a = origin[0]
        b = origin[1]
        c = origin[2]
        u = axis[0]
        v = axis[1]
        w = axis[2]
        # Set some intermediate values.
        u2 = u * u
        v2 = v * v
        w2 = w * w
        cos_t = cos(theta)
        sin_t = sin(theta)
        l2 = u2 + v2 + w2
        lsq = sqrt(l2)

        # Build the matrix entries element by element.
        m11 = (u2 + (v2 + w2) * cos_t) / l2
        m12 = (u * v * (1 - cos_t) - w * lsq * sin_t) / l2
        m13 = (u * w * (1 - cos_t) + v * lsq * sin_t) / l2
        m14 = (
            a * (v2 + w2) - u * (b * v + c * w) + (u * (b * v + c * w) - a * (v2 + w2)) * cos_t + (b * w - c * v) * lsq * sin_t
        ) / l2

        m21 = (u * v * (1 - cos_t) + w * lsq * sin_t) / l2
        m22 = (v2 + (u2 + w2) * cos_t) / l2
        m23 = (v * w * (1 - cos_t) - u * lsq * sin_t) / l2
        m24 = (
            b * (u2 + w2) - v * (a * u + c * w) + (v * (a * u + c * w) - b * (u2 + w2)) * cos_t + (c * u - a * w) * lsq * sin_t
        ) / l2

        m31 = (u * w * (1 - cos_t) - v * lsq * sin_t) / l2
        m32 = (v * w * (1 - cos_t) + u * lsq * sin_t) / l2
        m33 = (w2 + (u2 + v2) * cos_t) / l2
        m34 = (
            c * (u2 + v2) - w * (a * u + b * v) + (w * (a * u + b * v) - c * (u2 + v2)) * cos_t + (a * v - b * u) * lsq * sin_t
        ) / l2

        return SymmOp(
            [
                [m11, m12, m13, m14],
                [m21, m22, m23, m24],
                [m31, m32, m33, m34],
                [0, 0, 0, 1],
            ]
        )

    @staticmethod
    def reflection(
        normal: Union[Tuple[float, float, float], List[float], np.ndarray],
        origin: Union[Tuple[float, float, float], List[float], np.ndarray] = (0, 0, 0),
    ):
        """
        Returns reflection symmetry operation.
        Args:
            normal (3x1 array): Vector of the normal to the plane of
                reflection.
            origin (3x1 array): A point in which the mirror plane passes
                through.
        Returns:
            SymmOp for the reflection about the plane
        """
        # Normalize the normal vector first.
        n = np.array(normal, dtype=float) / np.linalg.norm(normal)

        u, v, w = n

        translation = np.eye(4)
        translation[0:3, 3] = -np.array(origin)

        xx = 1 - 2 * u**2
        yy = 1 - 2 * v**2
        zz = 1 - 2 * w**2
        xy = -2 * u * v
        xz = -2 * u * w
        yz = -2 * v * w
        mirror_mat = [[xx, xy, xz, 0], [xy, yy, yz, 0], [xz, yz, zz, 0], [0, 0, 0, 1]]

        if np.linalg.norm(origin) > 1e-6:
            mirror_mat = np.dot(np.linalg.inv(translation), np.dot(mirror_mat, translation))
        return SymmOp(mirror_mat)

    @staticmethod
    def inversion(
        origin: Union[Tuple[float, float, float], List[float], np.ndarray] = (0, 0, 0),
    ):
        """
        Inversion symmetry operation about axis.
        Args:
            origin (3x1 array): Origin of the inversion operation. Defaults
                to [0, 0, 0].
        Returns:
            SymmOp representing an inversion operation about the origin.
        """
        mat = -np.eye(4)
        mat[3, 3] = 1
        mat[0:3, 3] = 2 * np.array(origin)
        return SymmOp(mat)

    @staticmethod
    def roto_reflection(
        axis: Union[Tuple[float, float, float], List[float], np.ndarray],
        angle: float,
        origin: Union[Tuple[float, float, float], List[float], np.ndarray] = (0, 0, 0),
    ):
        """
        Returns a roto-reflection symmetry operation
        Args:
            axis (3x1 array): Axis of rotation / mirror normal
            angle (float): Angle in degrees
            origin (3x1 array): Point left invariant by roto-reflection.
                Defaults to (0, 0, 0).
        Return:
            Roto-reflection operation
        """
        rot = SymmOp.from_origin_axis_angle(origin, axis, angle)
        refl = SymmOp.reflection(axis, origin)
        m = np.dot(rot.affine_matrix, refl.affine_matrix)
        return SymmOp(m)

    def as_xyz_string(self) -> str:
        """
        Returns a string of the form 'x, y, z', '-x, -y, z',
        '-y+1/2, x+1/2, z+1/2', etc. Only works for integer rotation matrices
        """
        # test for invalid rotation matrix
        if not np.all(np.isclose(self.rotation_matrix, np.round(self.rotation_matrix))):
            warnings.warn('Rotation matrix should be integer')

        return transformation_to_string(self.rotation_matrix, translation_vec=self.translation_vector, delim=', ')

    @staticmethod
    def from_xyz_string(xyz_string: str):
        """
        Args:
            xyz_string: string of the form 'x, y, z', '-x, -y, z',
                '-2y+1/2, 3x+1/2, z-y+1/2', etc.
        Returns:
            SymmOp
        """
        rot_matrix = np.zeros((3, 3))
        trans = np.zeros(3)
        toks = xyz_string.strip().replace(' ', '').lower().split(',')
        re_rot = re.compile(r'([+-]?)([\d\.]*)/?([\d\.]*)([x-z])')
        re_trans = re.compile(r'([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])')
        for i, tok in enumerate(toks):
            # build the rotation matrix
            for m in re_rot.finditer(tok):
                factor = -1 if m.group(1) == '-' else 1
                if m.group(2) != '':
                    factor *= float(m.group(2)) / float(m.group(3)) if m.group(3) != '' else float(m.group(2))
                j = ord(m.group(4)) - 120
                rot_matrix[i, j] = factor
            # build the translation vector
            for m in re_trans.finditer(tok):
                factor = -1 if m.group(1) == '-' else 1
                num = float(m.group(2)) / float(m.group(3)) if m.group(3) != '' else float(m.group(2))
                trans[i] = num * factor
        return SymmOp.from_rotation_and_translation(rot_matrix, trans)
