import numpy as np
from itertools import chain, combinations

"""""choquet_integral takes a vector f which represent the fonction to integrate and a vector mu which represent
the fuzzy measure of each element of f. The  function computes the Choquet intergal"""


def choquet_integral(f, mu):
    n = len(mu)
    mu_sort = np.sort(mu)[::-1]
    integral = 0
    for i in range(n):
        integral += mu_sort[i] * np.sum(np.minimum(f, mu_sort[i]))
        f = np.maximum(f - mu_sort[i], 0)
    return integral


""""Examples of RIM quantifiers which represent the quantifiers “more than 100 ∗ 𝑘%” and “at least 100 ∗ 𝑘%”"""


def more_than(k, p):
    if p > k:
        return 1
    else:
        return 0


def at_least(k, p):
    if p >= k:
        return 1
    else:
        return 0


""""Zadeh’s S-function and the quantifiers "most" and "some"""


def zadeh_function(a, b, p):
    if p <= a:
        return 0
    elif a <= p <= (a + b) / 2:
        return (2 * (p - a) ** 2) / (b - a) ** 2
    elif (a + b) / 2 <= p <= b:
        return 1 - ((2 * (p - b) ** 2) / (b - a) ** 2)
    else:
        return 1


def most(p):
    return zadeh_function(0.3, 0.9, p)


def some(p):
    return zadeh_function(0.1, 0.4, p)


"""""The sets A_min and A_max of Definition 2.17"""


def A_min(gamma, val):
    if gamma > 0:
        if val > 0.5 * (1 + gamma):
            return 1
        elif 0.5 * (1 - gamma) < val < 0.5 * (1 + gamma):
            return 0.5
        else:
            return 0
    else:
        if val > 0.5:
            return 1
        elif val == 0.5:
            return 0.5
        else:
            return 0


def A_max(gamma, val):
    if gamma > 0:
        if val > 0.5 * (1 - gamma):
            return 0.5
        else:
            return 0
    else:
        if val >= 0.5:
            return 0.5
        else:
            return 0


""""The three-valued cut of 𝐴 at 𝛾"""


def three_valued_cut(gamma, val):
    return max(A_min(gamma, val), A_max(gamma, val))


""""The generalized fuzzy median"""


def gen_fuzzy_median(a, b):
    if min(a, b) > 0.5:
        return min(a, b)
    elif max(a, b) < 0.5:
        return max(a, b)
    else:
        return 0.5


""""The following functions are for the computation of the top and the bottom w.r.t gamma.
-Start by defining your own semi-fuzzy quantifier.
-The parameters of this function are the value of gamma, the set X which is a list(or numpy table), 
the set of fuzzy sets which is a list which contains a list of values of each fuzzy sets"""


def quantifier(sets):
    return 0


def powerset(s):
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def complement(x, y):
    return [el for el in y if el not in x]


def combination_tables(tab):
    n = len(tab)
    combinaisons = []

    def gen_com(i, combinaison):
        if i == n:
            combinaisons.append(combinaison)
            return
        for element in tab[i]:
            new_combinaison = combinaison.copy()
            new_combinaison.append(element)
            gen_com(i + 1, new_combinaison)

    gen_com(0, [])
    return combinaisons


def top_bottom(gamma, set, fuzzy_sets):
    tab_A_min = []
    tab_A_max = []
    tab_B = []
    for i in range(len(set)):
        tab_A_min.append([])
        tab_A_max.append([])

        if gamma > 0:
            for j in range(len(set)):
                if fuzzy_sets[i][j] >= 0.5 * (gamma + 1):
                    tab_A_min[i].append(set[j])
                if fuzzy_sets[i][j] > 0.5:
                    tab_A_max[i].append(set[j])
        if gamma == 0:
            for j in range(len(set)):
                if fuzzy_sets[i][j] > 0.5:
                    tab_A_min[i].append(set[j])
                if fuzzy_sets[i][j] >= 0.5:
                    tab_A_max[i].append(set[j])
    complement_tab = []
    powerset_table = []
    for i in range(len(tab_A_max)):
        complement_tab.append(complement(tab_A_min[i], tab_A_max[i]))
    for i in range(len(complement_tab)):
        powerset_table.append(list(powerset(complement(tab_A_min[i], tab_A_max[i]))))
    tab_B_final = []
    for i in range(len(tab_A_min)):
        tab_B.append([])
        for element in tab_A_min[i]:
            tab_B[i].append(element)

    tab = []
    for i in range(len(powerset_table)):
        tab.append([])
        for j in powerset_table[i]:
            tab[i].append(j)
    tab_B_save = tab_B
    for i in range(len(tab_B)):
        tab_B_final.append([])
        for k in range(len(tab[i])):
            for j in tab[i][k]:
                tab_B[i].append(j)
            tab_B_final[i].append(tab_B[i])
    for i in range(len(tab_B_save)):
        tab_B_final[i].append(tab_B_save[i])
    tab_combination = combination_tables(tab_B_final)
    val = []
    for i in range(len(tab_combination)):
        val.append(quantifier(tab_combination[i]))
    return max(val), min(val)


def top(gamma, set, fuzzy_sets):
    return top_bottom(gamma, set, fuzzy_sets)[0]


def bottom(gamma, set, fuzzy_sets):
    return top_bottom(gamma, set, fuzzy_sets)[1]


""""The quantification fuzzy mechanism:
-start by defining your own semi-fuzzy quantifier
-The parameters of this functions are as the parameters of the function top or bottom
"""


def QFM(gamma, set, fuzzy_sets):
    if bottom(gamma, set, fuzzy_sets) > 0.5:
        return bottom(gamma, set, fuzzy_sets)
    elif top(gamma, set, fuzzy_sets) < 0.5:
        return top(gamma, set, fuzzy_sets)
    else:
        return 0.5


"""""The following function compute the integral of Definition 2.21 using the trapezium method
-Define your own quantifier which will be used in the top and the buttom function
-set the parameter n which represents the number of subdividons of the interval [0,1]
-The parameters set and fuzzy_sets are defined as previously """


def QFM_OWA(set, fuzzy_sets, n):
    h = 1 / n
    result = 0.5 * ((top(0, set, fuzzy_sets)+bottom(0, set, fuzzy_sets))/2 + (top(1, set, fuzzy_sets)+bottom(1, set, fuzzy_sets))/2)
    for i in range(1, n):
        result += (top(i*h, set, fuzzy_sets)+bottom(i*h, set, fuzzy_sets))/2
    result *= h
    return result