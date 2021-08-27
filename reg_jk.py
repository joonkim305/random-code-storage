from sympy import symbols, solve
import numpy as np

def linear_reg(df_x, df_y, d):

    df_x = df_x.astype('float64')
    df_y = df_y.astype('float64')
    arr = [[0 for x in range(d + 1)] for y in range(d + 1)]

    for i in range(0, d + 1):
        for j in range(0, d + 1):
            arr[i][j] = sum(df_x ** (i + j))

    arr_inv = np.linalg.inv(arr)
    y_arr = [0 for x in range(d + 1)]
    for i in range(0, d + 1):
        y_arr[i] = sum(df_x ** i * df_y)

    coef = [0 for x in range(d + 1)]

    for i in range(0, d + 1):
        coef[i] = sum(y_arr * arr_inv[:][i])

    return coef

def solve_dx(coef):

    x = symbols('x', real=True)
    expr = 0
    coef_d = []
    for i in range(1, len(coef)):
        expr += i*coef[i] * (x ** (i-1))
        coef_d.append(i*coef[i])
    sol = solve(expr)
    return coef_d, sol

def local_min_max(coef, x_start, x_end):
    x_o = symbols('x_o', real=True)

    expr_org = 0
    for i in range(0, len(coef)):
        expr_org += coef[i] * x_o ** i

    coef_d, sol = solve_dx(coef)
    coef_d2, sol2 = solve_dx(coef_d)

    x = symbols('x', real=True)

    expr = 0
    for i in range(0, len(coef_d2)):
        expr += coef_d2[i] * x ** i

    local_min = []
    local_max = []
    abs_min = x_start

    abs_min_value = expr_org.subs(x_o, x_start)
    abs_max = x_start
    abs_max_value = expr_org.subs(x_o, x_start)
    for s in sol:

        if (s >= x_start) and (s <= x_end):
            d2_sol = expr.subs(x, s)

            if d2_sol > 0:
                local_min.append(s)
                if expr_org.subs(x_o, s) < abs_min_value:
                    abs_min_value = expr_org.subs(x_o, s)
                    abs_min = s
            else:
                local_max.append(s)
                if expr_org.subs(x_o, s) > abs_max_value:
                    abs_max_value = expr_org.subs(x_o, s)
                    abs_max = s

    return local_min, local_max, abs_min, abs_max