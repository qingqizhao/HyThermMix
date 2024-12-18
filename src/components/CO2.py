# from unittest import TestCase

# from scipy import arccosh, exp
from numpy import arccosh, exp


# from lib import unidades
# from lib.meos import MEoS


class CO2(object):
    """Multiparameter equation of state for carbon dioxide"""
    name = "carbon dioxide"
    CASNumber = "124-38-9"
    formula = "CO2"
    synonym = "R-744"
    _refPropName = "CO2"
    _coolPropName = "CarbonDioxide"
    # rhoc = unidades.Density(467.6)
    A = -2.63706E-01
    B = -8.82637E-01
    C = -1.19515E-02
    D = 2.21305E-03
    rhocL = 10.624978698
    rhoc = 467.6
    # Tc = unidades.Temperature(304.1282)
    Tc = 304.1282
    # Pc = unidades.Pressure(7.3773, "MPa")
    Pc = 7377.3 * 1e3  # "kPa"
    M = 44.0098  # g/mol
    # Tt = unidades.Temperature(216.592)
    Tt = 216.592
    # Tb = unidades.Temperature(194.686)
    Tb = 194.686
    f_acent = 0.22394
    # momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 49

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1],
           "ao_pow": [8.37304456, -3.70454304],
           "ao_exp": [1.99427042, .62105248, .41195293, 1.04028922, .08327678],
           "titao": [3.15163, 6.11190, 6.77708, 11.32384, 27.08792]}

    Fi2 = {"ao_log": [1, 2.50002],
           "pow": [0, 1], "ao_pow": [11.925152758, -16.118762264],
           "ao_exp": [], "titao": [],
           "ao_sinh": [2.04452, 2.03366], "sinh": [3.022758166, 1.589964364],
           "ao_cosh": [-1.06044, 0.01393], "cosh": [-2.844425476, 1.12159609]}

    CP3 = {"ao": 3.5,
           "ao_exp": [2, 1, 1], "exp": [960.11, 1932, 3380.2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Kunz "
                    "and Wagner (2004)",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 1100., "Pmax": 800000.0, "rhomax": 37.24,

        "nr1": [0.52646564804653, -0.14995725042592e1, 0.27329786733782,
                0.12949500022786],
        "d1": [1, 1, 2, 3],
        "t1": [0, 1.25, 1.625, 0.375],

        "nr2": [0.15404088341841, -0.58186950946814, -0.18022494838296,
                -0.095389904072812, -0.80486819317679e-2, -0.03554775127309,
                -0.28079014882405, -0.82435890081677e-1, 0.10832427979006e-1,
                -0.67073993161097e-2, -0.46827907600524e-2, -0.028359911832177,
                0.19500174744098e-1, -0.21609137507166, 0.43772794926972,
                -0.22130790113593, 0.15190189957331e-1, -0.15380948953300e-1],
        "d2": [3, 3, 4, 5, 6, 6, 1, 4, 1, 1, 3, 3, 4, 5, 5, 5, 5, 5],
        "t2": [0.375, 1.375, 1.125, 1.375, 0.125, 1.625, 3.75, 3.5, 7.5, 8, 6,
               16, 11, 24, 26, 28, 24, 26],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 6, 6],
        "gamma2": [1] * 18}

    eq = GERG