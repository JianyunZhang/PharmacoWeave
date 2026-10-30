"""Small standard-library linear-algebra and Student-t helpers.

These routines keep nonlinear least-squares identifiability checks and
asymptotic confidence intervals independent of heavy external numerical stacks.
"""

import math

_T_TWO_SIDED_95: dict[int, float] = {
    1: 12.7062047364,
    2: 4.3026527299,
    3: 3.1824463053,
    4: 2.7764451052,
    5: 2.5705818356,
    6: 2.4469118511,
    7: 2.3646242516,
    8: 2.3060041352,
    9: 2.2621571628,
    10: 2.2281388520,
    11: 2.2009851601,
    12: 2.1788128297,
    13: 2.1603686565,
    14: 2.1447866879,
    15: 2.1314495456,
    16: 2.1199052992,
    17: 2.1098155778,
    18: 2.1009220402,
    19: 2.0930240544,
    20: 2.0859634473,
    21: 2.0796138447,
    22: 2.0738730679,
    23: 2.0686576104,
    24: 2.0638985616,
    25: 2.0595385528,
    26: 2.0555294386,
    27: 2.0518305165,
    28: 2.0484071418,
    29: 2.0452296421,
    30: 2.0422724563,
}


def student_t_two_sided_quantile(df: int, ci_level: float = 0.95) -> float:
    """Return the upper critical value t_{(1 + ci_level)/2, df}."""

    if df < 1:
        raise ValueError("Degrees of freedom must be at least 1.")
    if not 0.0 < ci_level < 1.0:
        raise ValueError("ci_level must be strictly between 0 and 1.")
    if abs(ci_level - 0.95) <= 1e-12 and df in _T_TWO_SIDED_95:
        return _T_TWO_SIDED_95[df]
    p = 0.5 + 0.5 * ci_level
    if df == 1:
        return math.tan(math.pi * (p - 0.5))
    if df == 2:
        return (2.0 * p - 1.0) / math.sqrt(2.0 * p * (1.0 - p))
    z = _normal_quantile(p)
    z2 = z * z
    z3 = z2 * z
    z5 = z3 * z2
    z7 = z5 * z2
    z9 = z7 * z2
    nu = float(df)
    g1 = (z3 + z) / 4.0
    g2 = (5.0 * z5 + 16.0 * z3 + 3.0 * z) / 96.0
    g3 = (3.0 * z7 + 19.0 * z5 + 17.0 * z3 - 15.0 * z) / 384.0
    g4 = (79.0 * z9 + 776.0 * z7 + 1482.0 * z5 - 1920.0 * z3 - 945.0 * z) / 92160.0
    return z + g1 / nu + g2 / (nu * nu) + g3 / (nu**3) + g4 / (nu**4)


def solve_linear_system(matrix: list[list[float]], rhs: list[float]) -> list[float] | None:
    """Solve A x = b by Gaussian elimination with partial pivoting."""

    n = len(rhs)
    if n == 0 or len(matrix) != n or any(len(row) != n for row in matrix):
        return None
    aug = [[*map(float, row), float(b)] for row, b in zip(matrix, rhs, strict=True)]
    if any(not math.isfinite(val) for row in aug for val in row):
        return None
    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        pivot = aug[pivot_row][col]
        row_scale = max(1.0, max(abs(val) for val in aug[pivot_row][:n]))
        if abs(pivot) <= 1e-28 * row_scale:
            return None
        if pivot_row != col:
            aug[col], aug[pivot_row] = aug[pivot_row], aug[col]
        for r in range(col + 1, n):
            factor = aug[r][col] / aug[col][col]
            aug[r][col] = 0.0
            for c in range(col + 1, n + 1):
                aug[r][c] -= factor * aug[col][c]
    solution = [0.0] * n
    for r in range(n - 1, -1, -1):
        total = aug[r][n] - sum(aug[r][c] * solution[c] for c in range(r + 1, n))
        value = total / aug[r][r]
        if not math.isfinite(value):
            return None
        solution[r] = value
    return solution


def invert_symmetric(matrix: list[list[float]]) -> list[list[float]] | None:
    """Invert a small symmetric positive-definite matrix, or return None."""

    n = len(matrix)
    if n == 0:
        return None
    columns: list[list[float]] = []
    for index in range(n):
        unit = [1.0 if i == index else 0.0 for i in range(n)]
        col = solve_linear_system(matrix, unit)
        if col is None:
            return None
        columns.append(col)
    inv = [[0.5 * (columns[j][i] + columns[i][j]) for j in range(n)] for i in range(n)]
    if any(inv[i][i] <= 0.0 or not math.isfinite(inv[i][i]) for i in range(n)):
        return None
    return inv


def normalized_curvature_ratio(hessian: list[list[float]]) -> float:
    """Return min/max eigenvalue of the correlation-scaled Gram matrix."""

    n = len(hessian)
    if n == 0 or any(len(row) != n for row in hessian):
        return 0.0
    diag = [hessian[i][i] for i in range(n)]
    if any(not math.isfinite(val) or val <= 1e-24 for val in diag):
        return 0.0
    scales = [math.sqrt(val) for val in diag]
    work = [
        [1.0 if i == j else hessian[i][j] / (scales[i] * scales[j]) for j in range(n)]
        for i in range(n)
    ]
    if any(not math.isfinite(val) for row in work for val in row):
        return 0.0
    eigenvalues = _jacobi_eigenvalues(work)
    largest = max(eigenvalues)
    smallest = min(eigenvalues)
    if largest <= 1e-24 or smallest <= 0.0:
        return 0.0
    return smallest / largest


def _jacobi_eigenvalues(symmetric: list[list[float]]) -> list[float]:
    n = len(symmetric)
    a = [row[:] for row in symmetric]
    for _sweep in range(25):
        max_off = 0.0
        p, q = 0, 1
        for i in range(n):
            for j in range(i + 1, n):
                if abs(a[i][j]) > max_off:
                    max_off = abs(a[i][j])
                    p, q = i, j
        if max_off <= 1e-15:
            break
        diff = a[q][q] - a[p][p]
        if abs(diff) <= 1e-30:
            t = 1.0
        else:
            tau = diff / (2.0 * a[p][q])
            t = math.copysign(1.0 / (abs(tau) + math.sqrt(1.0 + tau * tau)), tau)
        c = 1.0 / math.sqrt(1.0 + t * t)
        s = t * c
        app = a[p][p] - t * a[p][q]
        aqq = a[q][q] + t * a[p][q]
        a[p][p] = app
        a[q][q] = aqq
        a[p][q] = a[q][p] = 0.0
        for k in range(n):
            if k in {p, q}:
                continue
            apk = a[p][k]
            aqk = a[q][k]
            a[p][k] = a[k][p] = c * apk - s * aqk
            a[q][k] = a[k][q] = s * apk + c * aqk
    return [a[i][i] for i in range(n)]


def _normal_quantile(p: float) -> float:
    a = (
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    )
    b = (
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    )
    c = (
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    )
    d = (
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    )
    p_low = 0.02425
    p_high = 1.0 - p_low
    if p < p_low:
        q = math.sqrt(-2.0 * math.log(p))
        num = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q) + c[5]
        den = ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q) + 1.0
        return num / den
    if p <= p_high:
        q = p - 0.5
        r = q * q
        num = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
        den = (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r) + 1.0
        return num / den
    q = math.sqrt(-2.0 * math.log(1.0 - p))
    num = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q) + c[5]
    den = ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q) + 1.0
    return -num / den
