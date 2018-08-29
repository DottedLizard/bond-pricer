class Solver(object):
    """Class provides root-finding methods for f(x)"""

    @staticmethod
    def bisection(func, a, b, target_tol, interval_tol, n_max):

        a = float(a)
        b = float(b)

        Solver.__validate_endpoints(a, b)

        if not Solver.__validate_axis_intersection(func, a, b):
            return None

        n = 0

        while n < n_max:

            c = (a + b) / 2
            c_val = func(c)

            if Solver.__within_tol(c_val, target_tol) or (((b - a) / 2) < interval_tol):
                return c

            n += 1

            if Solver.__sign(c_val) == Solver.__sign(func(a)):
                a = c
            else:
                b = c

        return None

    @staticmethod
    def __validate_endpoints(a, b):

        if a >= b:
            raise ValueError("a must be less than b")

    @staticmethod
    def __validate_axis_intersection(func, a, b):
        return not (((func(a) >= 0) and (func(b) >= 0)) or ((func(a) <= 0) and (func(b) <= 0)))

    @staticmethod
    def __within_tol(val, tol):
        return (val < tol) and (val > -tol)

    @staticmethod
    def __sign(x):

        if x > 0:
            return 1
        elif x < 0:
            return -1
        elif x == 0:
            return 0
        else:
            return x
