from math import exp


# from unittest import TestCase

# from lib import unidades
# from lib.meos import MEoS


class C2(object):
    """Multiparameter equation of state for ethane"""
    name = "ethane"
    CASNumber = "74-84-0"
    formula = "CH3CH3"
    synonym = "R-170"
    _refPropName = "ETHANE"
    _coolPropName = "Ethane"
    # rhoc = unidades.Density(206.18)
    A = -1.87015E-01
    B = -9.98938E-01
    C = -4.45252E-02
    D = -3.07766E-03
    rhocL = 6.870854540
    rhoc = 206.18
    # Tc = unidades.Temperature(305.322)
    Tc = 305.322
    # Pc = unidades.Pressure(4872.2, "kPa")
    Pc = 4872.2 * 1e3  # "kPa"
    M = 30.06904  # g/mol
    # Tt = unidades.Temperature(90.368)
    Tt = 90.368
    # Tb = unidades.Temperature(184.569)
    Tb = 184.569
    f_acent = 0.0995
    # momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 3
    # _Tr = unidades.Temperature(295.159630)
    _Tr = 295.159630
    # _rhor = unidades.Density(207.557649)
    _rhor = 207.557649
    _w = 0.095234716

    Fi1 = {"R": 8.314472,
           "ao_log": [1, 3.003039265],
           "pow": [0, 1],
           "ao_pow": [9.212802589, -4.68224855],
           "ao_exp": [1.117433359, 3.467773215, 6.941944640, 5.970850948],
           "titao": [1.4091052332, 4.0099170712, 6.5967098342, 13.9798102659]}

    Fi2 = {"R": 8.31451,
           "ao_log": [1, 3.00263],
           "pow": [0, 1],
           "ao_pow": [24.675437527, -77.42531376],
           "ao_exp": [], "titao": [],
           "ao_sinh": [4.33939, 13.1974],
           "sinh": [559.314 / Tc, 1031.38 / Tc],
           "ao_cosh": [1.23722, -6.01989],
           "cosh": [223.284 / Tc, 1071.29 / Tc]}

    Fi3 = {"ao_log": [1, 3.8159476],
           "pow": [0, -1. / 3, -2. / 3, -1],
           "ao_pow": [-23.446765, 8.6021299, -3.3075735, -0.55956678],
           "ao_exp": [5.0722267], "titao": [5.5074874]}

    CP5 = {"ao": 9.9507922459,
           "an": [-6.9341406909e5, 3.1534834135e4, -6.103375287e2,
                  -2.8657877948e-2, 9.0922897821e-5, -5.2750109915e-8],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [-1.4243593411e1], "exp": [3000]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethane of Kunz and Wagner"
                    " (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 675.0, "Pmax": 900000.0, "rhomax": 22.419,

        "nr1": [0.63596780450714, -0.17377981785459e1, 0.28914060926272,
                -0.33714276845694, 0.22405964699561e-1, 0.15715424886913e-1],
        "d1": [1, 1, 2, 2, 4, 4],
        "t1": [0.125, 1.125, 0.375, 1.125, 0.625, 1.5],

        "nr2": [0.11450634253745, 0.10612049379745e1, -0.12855224439423e1,
                0.39414630777652, 0.31390924682041, -0.21592277117247e-1,
                -0.21723666564905, -0.28999574439489, 0.42321173025732,
                0.46434100259260e-1, -0.13138398329741, 0.11492850364368e-1,
                -0.33387688429909e-1, 0.015183171583644, -0.47610805647657e-2,
                0.46917166277885e-1, -0.039401755804649, -0.32569956247611e-2],
        "d2": [1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7],
        "t2": [0.625, 2.625, 2.75, 2.125, 2, 1.75, 4.5, 4.75, 5, 4, 4.5, 7.5,
               14, 11.5, 26, 28, 30, 16],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6],
        "gamma2": [1] * 18}

    eq = GERG