from math import exp, log, log10


# from unittest import TestCase

# from lib import unidades
# from lib.meos import MEoS


class H2(object):
    """Multiparameter equation of state for hydrogen (normal)"""
    name = "hydrogen"
    CASNumber = "1333-74-0"
    formula = "H2"
    synonym = "R-702"
    _refPropName = "HYDROGEN"
    _coolPropName = "Hydrogen"
    # rhoc = unidades.Density(31.26226704)
    A = 1.25418E-01
    B = -5.08704E-01
    C = 8.12269E-01
    D = 2.04251E-02
    rhocL = 14.94
    rhoc = 30.1172472
    # Tc = unidades.Temperature(33.145)
    Tc = 33.19
    # Pc = unidades.Pressure(1296.4, "kPa")
    Pc = 1296.4 * 1e3  # "kPa"
    M = 2.01588  # g/mol
    #     Tt = unidades.Temperature(13.957)
    #     Tb = unidades.Temperature(20.369)
    Tt = 13.957
    Tb = 20.369
    f_acent = -0.219
    # momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 1

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-1.4579856475, 1.888076782],
           "ao_exp": [1.616, -0.4117, -0.792, 0.758, 1.217],
           "titao": [16.0205159149, 22.6580178006, 60.0090511389,
                     74.9434303817, 206.9392065168]}

    Fi2 = {"ao_log": [1, 1.47906],
           "pow": [0, 1],
           "ao_pow": [13.796443393, -175.864487294],
           "ao_sinh": [0.95806, 1.56039], "sinh": [6.891654113, 49.76529075],
           "ao_cosh": [0.45444, -1.3756], "cosh": [9.84763483, 50.367279301]}

    CP1 = {"ao": 0.72480209e3,
           "an": [0.12155215e11, -0.36396763e10, 0.43375265e9, -0.23085817e8,
                  -0.38680927e4, 0.88240136e5, -0.78587085e4, -0.18426806e3,
                  0.21801550e2, -0.13051820e1, 0.21003175e-1, 0.23911604e-2,
                  -0.18240547e-3, 0.56149561e-5, -0.73803310e-7,
                  0.66357755e-11],
           "pow": [-7, -6, -5, -4, -3, -2, -1.001, 0.5, 1, 1.5, 2, 2.5, 3, 3.5,
                   4, 5]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 121000.0, "rhomax": 38.148,

        "nr1": [0.53579928451252e1, -0.62050252530595e1, 0.13830241327086,
                -0.71397954896129e-1, 0.15474053959733e-1],
        "d1": [1, 1, 2, 2, 4],
        "t1": [0.5, 0.625, 0.375, 0.625, 1.125],

        "nr2": [-0.14976806405771, -0.26368723988451e-1, 0.56681303156066e-1,
                -0.60063958030436e-1, -0.45043942027132, 0.42478840244500,
                -0.021997640827139, -0.01049952137453, -0.28955902866816e-2],
        "d2": [1, 5, 5, 5, 1, 1, 2, 5, 1],
        "t2": [2.625, 0, 0.25, 1.375, 4, 4.25, 5, 8, 8],
        "c2": [1, 1, 1, 1, 2, 2, 3, 3, 5],
        "gamma2": [1] * 9}

    eq = GERG