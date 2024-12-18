from math import exp, pi

from scipy.constants import Boltzmann


class Ar(object):
    """Multiparamter equation of state for argon"""
    name = "argon"
    CASNumber = "7440-37-1"
    formula = "Ar"
    synonym = "R-740"
    _refPropName = "ARGON"
    _coolPropName = "Argon"
    A = -6.55869E-02
    B = -9.61670E-01
    C = 4.85341E-02
    D = 7.28174E-03
    rhocL = 13.407429659
    rhoc = 535.5999973
    Tc = 150.687
    Pc = 4863 * 1e3  # "kPa"
    M = 39.948  # g/mol
    Tt = 83.8058
    Tb = 87.302
    f_acent = -0.00219
    # momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 98
    _Tr = 147.707801
    _rhor = 540.014968
    _w = 0.000305675

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [8.31666243, -4.94651164]}

    CP1 = {"ao": 2.5}

    Fi2 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [8.3166315, -4.9465026]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for argon of Kunz and Wagner "
                    "(2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": 83.8, "Tmax": 700., "Pmax": 1000000.0, "rhomax": 50.65,

        "nr1": [0.85095714803969, -0.24003222943480e1, 0.54127841476466,
                0.16919770692538e-1, 0.68825965019035e-1, 0.21428032815338e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.17429895321992, -0.033654495604194, -0.13526799857691,
                -0.016387350791552, -0.024987666851475, 0.0088769204815709],
        "c2": [1, 1, 2, 2, 3, 3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "gamma2": [1] * 6,

        "nr3": [],
        "nr4": []}

    # eq = tegeler, younglove, GERG, stewart, shortSpan