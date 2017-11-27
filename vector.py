"""
Vector class from Udacity Linear Algebra Refresher course
"""

import math
from decimal import Decimal as D
from decimal import getcontext

getcontext().prec = 30
EPSILON = 1e-8
DIGIT_KEPT = 8

class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([D(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(tuple([round(x, DIGIT_KEPT) for x in
                                          self.coordinates]))


    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def add(self, v2):
        if self.dimension == v2.dimension:
            return Vector([self.coordinates[i] + D(v2.coordinates[i]) for i in
                           range(v1.dimension)])
        return None

    def minus(self, v2):
        if self.dimension == v2.dimension:
            return Vector([self.coordinates[i] - D(v2.coordinates[i]) for i in
                           range(self.dimension)])
        return None

    def scalar_multiply(self, s):
        return Vector([self.coordinates[i] * D(s) for i in range(
            self.dimension)])

    def magnitude(self):
        l2_squred = D(sum([self.coordinates[i]**2 for i in range(self.dimension)]))
        return math.sqrt(l2_squred)

    def normalization(self):
        try:
            return self.scalar_multiply(1.0 / self.magnitude())
        except ZeroDivisionError:
            raise Exception("l2 norm is zero. Cannot normalize")

    def dot_product(self, v):
        if v.dimension == self.dimension:
            return sum([x*D(y) for x,y in zip(self.coordinates,
                                              v.coordinates)])
        return None

    def get_angle(self, v, in_degree=False):
        norm_1 = D(self.magnitude())
        norm_2 = D(v.magnitude())
        if norm_1 == 0.0 or norm_2 == 0.0:
            raise Exception("l2 norm is zero for one or both vectors. "
                            "Either orthogonal or no angle definition.")
        else:
            angle_in_radian = math.acos(
                round(self.dot_product(v)/(norm_1 * norm_2), DIGIT_KEPT))
            if in_degree:
                return math.degrees(angle_in_radian)
            return angle_in_radian


    def is_zero_vector(self):
        return all(abs(element) < EPSILON for element in self.coordinates)


    def is_parallel(self, v):
        return (self.is_zero_vector() or v.is_zero_vector() or
                self.get_angle(v) == math.pi or self.get_angle(v) == 0.0)


    def is_orthogonal(self, v):
        if self.is_zero_vector() or v.is_zero_vector():
            return True
        elif abs(self.dot_product(v)) < EPSILON:
            return True
        else:
            return False

    def get_parallel_projection(self, b):
        try:
            # alternatively
            # v_magnitude = self.magnitude()
            # v_parallel_magnitude = D(math.cos(self.get_angle(b)) * v_magnitude)
            unit_basis = b.normalization()
            v_parallel_magnitude = self.dot_product(unit_basis)
            return unit_basis.scalar_multiply(v_parallel_magnitude)
        except Exception as e:
            raise e


    def get_orthogonal_projection(self, b):
        try:
            # alternatively
            # v_magnitude = self.magnitude()
            # v_orthogonal_magnitude = D(math.sin(self.get_angle(b)) * v_magnitude)
            unit_orthogonal_to_b = \
                self.minus(self.get_parallel_projection(b)).normalization()
            v_orthogonal_magnitude = self.dot_product(unit_orthogonal_to_b)
            return unit_orthogonal_to_b.scalar_multiply(v_orthogonal_magnitude)
        except Exception as e:
            raise e

    def cross_product(self, v):
        if self.dimension != 3 or v.dimension != 3:
            raise Exception("Not 3 dimentions in calculating cross "
                            "product.")
        else:
            if self.is_zero_vector() or v.is_zero_vector() or \
                    self.is_parallel(v):
                return Vector([0.0, 0.0, 0.0])
            else:
                return Vector([self.coordinates[1] * v.coordinates[2] -
                               v.coordinates[1] * self.coordinates[2],
                               -(self.coordinates[0] * v.coordinates[2] -
                                v.coordinates[0] * self.coordinates[2]),
                               self.coordinates[0] * v.coordinates[1] -
                               v.coordinates[0] * self.coordinates[1]])


    def get_parallelogram_area(self, v):
        try:
            crs_prd = self.cross_product(v)
            return D(crs_prd.magnitude())
        except Exception as e:
            raise e


    def get_triangle_area(self, v):
        try:
            return self.get_parallelogram_area(v)/D(2.0)
        except Exception as e:
            raise e





v1 = Vector([8.218, -9.341])
v2 = Vector([-1.129, 2.111])
print v1.add(v2)

v3 = Vector([7.119, 8.215])
v4 = Vector([-8.223, 0.878])
print v3.minus(v4)

v5 = Vector([1.671, -1.012, -0.318])
print v5.scalar_multiply(7.41)


v6 = Vector([-0.221, 7.437])
v7 = Vector([8.813, -1.331, -6.247])
print v6.magnitude(), v7.magnitude()

v8 = Vector([5.581, -2.136])
v9 = Vector([1.996, 3.108, -4.554])
print v8.normalization(), v9.normalization()

v10 = Vector([7.887, 4.138])
v11 = Vector([-8.802, 6.776])
print v10.dot_product(v11)


v12 = Vector([-5.955, -4.904, -1.874])
v13 = Vector([-4.496, -8.755, 7.103])
print v12.dot_product(v13)


v14 = Vector([3.183, -7.627])
v15 = Vector([-2.668, 5.319])
print v14.get_angle(v15)

v16 = Vector([7.35, 0.221, 5.188])
v17 = Vector([2.751, 8.259, 3.985])
print v16.get_angle(v17, True)

v18 = Vector([-7.579, -7.88])
v19 = Vector([22.737, 23.64])
print v18.is_parallel(v19), v18.is_orthogonal(v19)
# print v18.dot_product(v19)
# print [v18.coordinates[i]/v19.coordinates[i] for i in range(v18.dimension)]

v20 = Vector([-2.029, 9.97, 4.172])
v21 = Vector([-9.231, -6.639, -7.245])
print v20.is_parallel(v21), v20.is_orthogonal(v21)

v22 = Vector([-2.328, -7.284, -1.214])
v23 = Vector([-1.821, 1.072, -2.94])
print v22.is_parallel(v23), v22.is_orthogonal(v23)

v24 = Vector([2.118, 4.827])
v25 = Vector([0, 0])
print v24.is_parallel(v25), v24.is_orthogonal(v25)

v26 = Vector([3.039, 1.879])
v27 = Vector([0.825, 2.036])
print v26.get_parallel_projection(v27)

v28 = Vector([-9.88, -3.264, -8.159])
v29 = Vector([-2.155, -9.353, -9.473])
print v28.get_orthogonal_projection(v29)

v30 = Vector([3.009, -6.172, 3.692, -2.51])
v31 = Vector([6.404, -9.144, 2.759, 8.718])

print v30.get_parallel_projection(v31),\
    v30.get_orthogonal_projection(v31)

v32 = Vector([8.462, 7.893, -8.187])
v33 = Vector([6.984, -5.975, 4.778])
print v32.cross_product(v33)

v34 = Vector([-8.987, -9.838, 5.031])
v35 = Vector([-4.268, -1.861, -8.866])
print v34.get_parallelogram_area(v35)

v36 = Vector([1.5, 9.547, 3.691])
v37 = Vector([-6.007, 0.124, 5.772])
print v36.get_triangle_area(v37)





