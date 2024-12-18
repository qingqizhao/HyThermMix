from math import exp



class H2S(object):
    """Multiparameter equation of state for hydrogen sulfide"""
    name = "hydrogen sulfide"
    CASNumber = "7783-06-4"
    formula = "H2S"
    synonym = ""
    _refPropName = "H2S"
    _coolPropName = "HydrogenSulfide"
    A = -5.87760E-02
    B = -8.44317E-01
    C = -1.35260E-02
    D = -3.72996E-03
    rhocL = 10.19
    rhoc = 347.2841672
    Tc = 373.1
    Pc = 9000.0*1e3 # "kPa"
    M = 34.08088  # g/mol
    Tt = 187.7
    Tb = 212.85
    f_acent = 0.1005
    #momentoDipolar = unidades.DipoleMoment(0.97, "Debye")
    id = 50

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1.5],
           "ao_pow": [-4.0740770957, 3.7632137341, -0.002753352822675789],
           "ao_exp": [1.1364, 1.9721],
           "titao": [1823/Tc, 3965/Tc]}

    Fi2 = {"ao_log": [1, 3],
           "pow": [0, 1], "ao_pow": [9.336197742, -16.266508995],
           "ao_sinh": [3.11942], "sinh": [4.914580541],
           "ao_cosh": [1.00243], "cosh": [2.27065398]}

    Fi3 = {"ao_log": [1, 3.],
           "pow": [0, 1], "ao_pow": [7.881037, -3.20986],
           "ao_exp": [0.9767422, 2.151898], "titao": [4.506266, 10.15526]}

    CP3 = {"ao": 4.1012105,
           "an": [-1.6720073e-3, 7.5303152e-6, -0.62421053e-8, 0.18098453e-11],
           "pow": [1, 2, 3, 4]}

    f = 8.3159524/4.184
    CP4 = {"ao": 7.9468*f,
           "ao_sinh": [-1.5769761e4*f, 1.3861204e7*f],
           "sinh": [4.33801e2, 1.48143e3],
           "ao_cosh": [2.0329947e6*f, -3.5044957e6*f],
           "cosh": [8.43792e2, 1.10223e3]}


    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Kunz and "
                    "Wagner (2008).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 760.0, "Pmax": 170000.0, "rhomax": 29.12,

        "nr1":  [0.87641, -2.0367, 0.21634, -0.050199, 0.066994, 0.19076e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.20227, -0.45348e-2, -0.2223, -0.034714, -.014885, .74154e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    #eq = lemmon, sakoda, polt, starling, GERG