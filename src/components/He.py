from math import exp, log


class He(object):
    """Multiparameter equation of state for helium"""
    name = "helium"
    CASNumber = "7440-59-7"
    formula = "He"
    synonym = "R-704"
    _refPropName = "HELIUM"
    _coolPropName = "Helium"
    A = 2.54306E-01
    B = -1.44743E-01
    C = 8.90411E-01
    D = 2.03250E-03
    rhocL = 17.399
    rhoc = 69.6412722
    Tc = 5.1953
    Pc = 227.61 * 1e3  # "kPa"
    M = 4.002602  # g/mol
    Tt = 2.1768
    Tb = 4.2226
    f_acent = -0.385
    # momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 212

    CP1 = {"ao": 2.5}

    Fi2 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [13.628409737, -143.470759602]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for helium of Kunz and "
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

        "Tmin": Tt, "Tmax": 1500.0, "Pmax": 100000.0, "rhomax": 88.73,

        "nr1": [-0.45579024006737, 0.12516390754925e1, -0.15438231650621e1,
                0.20467489707221e-1],
        "d1": [1, 1, 1, 4],
        "t1": [0, 0.125, 0.75, 1.],

        "nr2": [-0.34476212380781, -0.20858459512787e-1, 0.16227414711778e-1,
                -0.057471818200892, 0.19462416430715e-1, -0.33295680123020e-1,
                -0.10863577372367e-1, -0.22173365245954e-1],
        "d2": [1, 3, 5, 5, 5, 2, 1, 2],
        "t2": [0.75, 2.625, 0.125, 1.25, 2., 1., 4.5, 5.],
        "c2": [1, 1, 1, 1, 1, 2, 3, 3],
        "gamma2": [1] * 8,

        "nr3": [],
        "nr4": []}

    # eq = ortiz, mccarty, MBWR, GERG
