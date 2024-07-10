# CUBIC ROOT SOLVER

# Date Created   :    24.05.2017
# Created by     :    Shril Kumar [(shril.iitdhn@gmail.com),(github.com/shril)] &
#                     Devojoyti Halder [(devjyoti.itachi@gmail.com),(github.com/devojoyti)]

# Project        :    Classified 
# Use Case       :    Instead of using standard numpy.roots() method for finding roots,
#                     we have implemented our own algorithm which is ~10x faster than
#                     in-built method.

# Algorithm Link :    www.1728.org/cubic2.htm

# This script (Cubic Equation Solver) is an independent program for computation of roots of Cubic Polynomials. This script, however,
# has no relation with original project code or calculations. It is to be also made clear that no knowledge of it's original project 
# is included or used to device this script. This script is complete freeware developed by above signed users, and may further be
# used or modified for commercial or non-commercial purpose.


# Libraries imported for fast mathematical computations.
import math
import numpy as np

# Main Function takes in the coefficient of the Cubic Polynomial
# as parameters and it returns the roots in form of numpy array.
# Polynomial Structure -> ax^3 + bx^2 + cx + d = 0

def solve(a, b, c, d):
    """
    Solves the cubic equation ax^3 + bx^2 + cx + d = 0 and returns its roots.

    This function computes the roots of a cubic polynomial equation using a custom algorithm 
    which is significantly faster than the standard numpy.roots() method. It can also handle 
    special cases of linear and quadratic equations when the higher degree coefficients are zero.

    Parameters:
    a (float): Coefficient of x^3
    b (float): Coefficient of x^2
    c (float): Coefficient of x
    d (float): Constant term

    Returns:
    numpy.ndarray: An array containing the roots of the polynomial. The roots can be real or complex numbers.

    Notes:
    - If both a and b are zero, the function treats the equation as a linear equation and returns the single root.
    - If a is zero, the function treats the equation as a quadratic equation and returns its two roots.
    - For cubic equations, it uses different methods depending on the discriminant:
      - If f, g, and h are zero, all roots are real and equal.
      - If h is less than or equal to zero, all three roots are real.
      - If h is greater than zero, there is one real root and two complex conjugate roots.

    Example Usage:
    >>> solve(1, -6, 11, -6)
    array([3., 2., 1.])

    >>> solve(0, 1, -3, 2)
    array([2., 1.])

    >>> solve(0, 0, 1, -2)
    array([2.])

    """
    if (a == 0 and b == 0):                     # Case for handling Liner Equation
        return np.array([(-d * 1.0) / c])                 # Returning linear root as numpy array.

    elif (a == 0):                              # Case for handling Quadratic Equations

        D = c * c - 4.0 * b * d                       # Helper Temporary Variable
        if D >= 0:
            D = math.sqrt(D)
            x1 = (-c + D) / (2.0 * b)
            x2 = (-c - D) / (2.0 * b)
        else:
            D = math.sqrt(-D)
            x1 = (-c + D * 1j) / (2.0 * b)
            x2 = (-c - D * 1j) / (2.0 * b)
            
        return np.array([x1, x2])               # Returning Quadratic Roots as numpy array.

    f = findF(a, b, c)                          # Helper Temporary Variable
    g = findG(a, b, c, d)                       # Helper Temporary Variable
    h = findH(g, f)                             # Helper Temporary Variable

    if f == 0 and g == 0 and h == 0:            # All 3 Roots are Real and Equal
        if (d / a) >= 0:
            x = (d / (1.0 * a)) ** (1 / 3.0) * -1
        else:
            x = (-d / (1.0 * a)) ** (1 / 3.0)
        return np.array([x, x, x])              # Returning Equal Roots as numpy array.

    elif h <= 0:                                # All 3 roots are Real

        i = math.sqrt(((g ** 2.0) / 4.0) - h)   # Helper Temporary Variable
        j = i ** (1 / 3.0)                      # Helper Temporary Variable
        k = math.acos(-(g / (2 * i)))           # Helper Temporary Variable
        L = j * -1                              # Helper Temporary Variable
        M = math.cos(k / 3.0)                   # Helper Temporary Variable
        N = math.sqrt(3) * math.sin(k / 3.0)    # Helper Temporary Variable
        P = (b / (3.0 * a)) * -1                # Helper Temporary Variable

        x1 = 2 * j * math.cos(k / 3.0) - (b / (3.0 * a))
        x2 = L * (M + N) + P
        x3 = L * (M - N) + P

        return np.array([x1, x2, x3])           # Returning Real Roots as numpy array.

    elif h > 0:                                 # One Real Root and two Complex Roots
        R = -(g / 2.0) + math.sqrt(h)           # Helper Temporary Variable
        if R >= 0:
            S = R ** (1 / 3.0)                  # Helper Temporary Variable
        else:
            S = (-R) ** (1 / 3.0) * -1          # Helper Temporary Variable
        T = -(g / 2.0) - math.sqrt(h)
        if T >= 0:
            U = (T ** (1 / 3.0))                # Helper Temporary Variable
        else:
            U = ((-T) ** (1 / 3.0)) * -1        # Helper Temporary Variable

        x1 = (S + U) - (b / (3.0 * a))
        x2 = -(S + U) / 2 - (b / (3.0 * a)) + (S - U) * math.sqrt(3) * 0.5j
        x3 = -(S + U) / 2 - (b / (3.0 * a)) - (S - U) * math.sqrt(3) * 0.5j

        return np.array([x1, x2, x3])           # Returning One Real Root and two Complex Roots as numpy array.


# Helper function to return float value of f.
def findF(a, b, c):
    return ((3.0 * c / a) - ((b ** 2.0) / (a ** 2.0))) / 3.0


# Helper function to return float value of g.
def findG(a, b, c, d):
    return (((2.0 * (b ** 3.0)) / (a ** 3.0)) - ((9.0 * b * c) / (a **2.0)) + (27.0 * d / a)) /27.0


# Helper function to return float value of h.
def findH(g, f):
    return ((g ** 2.0) / 4.0 + (f ** 3.0) / 27.0)
