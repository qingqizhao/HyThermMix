from math import log, exp


class nC4(object):
    """Multiparameter equation of state for n-butane"""
    name = "n-butane"
    CASNumber = "106-97-8"
    formula = "CH3-(CH2)2-CH3"
    synonym = "R-600"
    _refPropName = "BUTANE"
    _coolPropName = "n-Butane"
    # rhoc = unidades.Density(228.)
    A = -2.05100E-01
    B = -8.31992E-01
    C = -1.48650E-02
    D = -5.46152E-04
    rhocL = 3.920016792
    rhoc = 228.
    # Tc = unidades.Temperature(425.125)
    Tc = 425.125
    # Pc = unidades.Pressure(3796.0, "kPa")
    Pc = 3796.0 * 1e3  # "kPa"
    M = 58.1222  # g/mol
    # Tt = unidades.Temperature(134.895)
    Tt = 134.895
    # Tb = unidades.Temperature(272.660)
    Tb = 272.660
    f_acent = 0.201
    # momentoDipolar = unidades.DipoleMoment(0.05, "Debye")
    id = 6
    # _Tr = unidades.Temperature(406.785141)
    _Tr = 406.785141
    # _rhor = unidades.Density(230.384826)
    _rhor = 230.384826
    _w = 0.194240287

    Fi1 = {"ao_log": [1, 3.24680487],
           "pow": [0, 1],
           "ao_pow": [12.54882924, -5.46976878],
           "ao_exp": [5.54913289, 11.4648996, 7.59987584, 9.66033239],
           "titao": [0.7748404445, 3.3406025522, 4.9705130961, 9.9755537783]}

    Fi2 = {"ao_log": [1, 3.33944],
           "pow": [0, 1],
           "ao_pow": [20.884143364, -91.638478026],
           "ao_sinh": [9.44893, 24.4618],
           "sinh": [468.27 / Tc, 1914.1 / Tc],
           "ao_cosh": [6.89406, 14.7824],
           "cosh": [183.636 / Tc, 903.185 / Tc]}

    Fi3 = {"ao_log": [1, 3.240207],
           "pow": [0, 1],
           "ao_pow": [-5.404217, 4.91136],
           "ao_exp": [5.513671, 7.388450, 10.250630, 11.061010],
           "titao": [327.55988 / Tc, 1319.06935 / Tc,
                     4138.63184 / Tc, 1864.36783 / Tc]}

    CP4 = {"ao": -1.3491511376e1,
           "an": [3.8802310194e5, -1.5444296890e5, 2.8455082239e3,
                  6.6142595353e-2, -2.4307965028e-5, 1.5044248429e-10],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [-8.3933423467], "exp": [3000]}

    CP6 = {"ao": 0.801601 / 8.3143 * 58.124,
           "an": [0.655936e-3 / 8.3143 * 58.124, 0.12277e-4 / 8.3143 * 58.124,
                  -0.165626e-7 / 8.3143 * 58.124, 0.67736e-11 / 8.3143 * 58.124],
           "pow": [1, 2, 3, 4]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Kunz and "
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

        "Tmin": Tt, "Tmax": 575., "Pmax": 69000.0, "rhomax": 13.2,

        "nr1": [0.10626277411455e1, -0.28620951828350e1, 0.88738233403777,
                -0.12570581155345, 0.10286308708106, 0.25358040602654e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.32325200233982, -0.037950761057432, -0.32534802014452,
                -0.079050969051011, -0.020636720547775, 0.57053809334750e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1] * 6}

    eq = GERG