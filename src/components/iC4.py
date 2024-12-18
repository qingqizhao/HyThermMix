from math import exp


class iC4(object):
    """Multiparameter equation of state for isobutane"""
    name = "isobutane"
    CASNumber = "75-28-5"
    formula = "CH(CH3)3"
    synonym = "R-600a"
    _refPropName = "ISOBUTAN"
    _coolPropName = "IsoButane"
    A = -1.70080E-01
    B = -8.31823E-01
    C = -1.10934E-02
    D = -1.09953E-03
    rhocL = 3.860142940
    rhoc = 224.36
    Tc = 407.817
    Pc = 3629.0 * 1e3  # "kPa"
    M = 58.1222  # g/mol
    Tt = 113.73
    Tb = 261.401
    f_acent = 0.184
    # momentoDipolar = unidades.DipoleMoment(0.132, "Debye")
    id = 5
    _Tr = 390.355535
    _rhor = 228.302484
    _w = 0.178714317

    Fi1 = {"ao_log": [1, 3.05956619],
           "pow": [0, 1],
           "ao_pow": [11.60865546, -5.29450411],
           "ao_exp": [4.94641014, 4.09475197, 15.6632824, 9.73918122],
           "titao": [0.9512779015, 2.3878958853, 4.3469042691, 10.3688586351]}

    Fi2 = {"ao_log": [1, 3.06714],
           "pow": [0, 1],
           "ao_pow": [20.413726078, -94.467620036],
           "ao_sinh": [8.97575, 25.1423], "sinh": [1.074673199, 4.671261865],
           "ao_cosh": [5.25156, 16.1388], "cosh": [0.485556021, 2.19158348]}

    Fi3 = {"ao_log": [1, 3.059347],
           "ao_pow": [-6.026745, -5.035251],
           "pow": [0, 1],
           "ao_exp": [4.940314, 4.090139, 9.739581, 15.68832],
           "titao": [0.9508183, 2.383449, 10.38655, 4.347095]}

    CP4 = {"ao": -1.7231723278e1,
           "an": [1.7027919006e7, -4.7269724737e5, 4.7301406581e3,
                  5.8491344291e-2, 8.9440351886e-6, -1.8274599197e-8],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [-1.9283021962e1], "exp": [3000]}

    CP5 = {"ao": 4.06714,
           "ao_sinh": [8.97575, 25.1423], "sinh": [438.27, 1905.02],
           "ao_cosh": [5.25156, 16.1388], "cosh": [198.018, 893.765]}

    CP6 = {"ao": 0.397893 / 8.3143 * 58.124,
           "an": [0.412501e-2 / 8.3143 * 58.124, -0.196195e-6 / 8.3143 * 58.124,
                  0.380185e-8 / 8.3143 * 58.124, -0.523950e-11 / 8.3143 * 58.124],
           "pow": [1, 2, 3, 4]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Kunz and "
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

        "Tmin": Tt, "Tmax": 575.0, "Pmax": 35000.0, "rhomax": 12.9,

        "nr1": [0.10429331589100e1, -0.28184272548892e1, 0.86176232397850,
                -0.10613619452487, 0.98615749302134e-1, 0.23948208682322e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.30330004856950, -0.041598156135099, -0.29991937470058,
                -0.080369342764109, -0.029761373251151, 0.013059630303140],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1] * 6}

    # eq = buecker, younglove, GERG, miyamoto, shortSpan, polt, sun