from math import exp
#from unittest import TestCase

#from lib import unidades
#from lib.meos import MEoS


class CH4(object):
    """Multiparameter equation of state for methane"""
    name = "methane"
    CASNumber = "74-82-8"
    formula = "CH4"
    synonym = "R-50"
    _refPropName = "METHANE"
    _coolPropName = "Methane"
    #rhoc = unidades.Density(162.66)
    A = -3.14160e-02
    B = -7.22738E-01
    C = 1.13859E-01
    D = 1.16905E-02
    rhocL = 10.139342719
    rhoc = 162.66
    #Tc = unidades.Temperature(190.564)
    Tc =  190.564
    #Pc = unidades.Pressure(4599.2, "kPa")
    Pc = 4599.2*1e3 # "kPa"
    M = 16.0428  # g/mol
    #Tt = unidades.Temperature(90.694)
    Tt = 90.694
    #Tb = unidades.Temperature(111.667)
    Tb = 111.667
    f_acent = 0.01142
    #momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 2
    #_Tr = unidades.Temperature(186.659809)
    _Tr = 186.659809
    #_rhor = unidades.Density(163.413536)
    _rhor = 163.413536
    _w = 0.010528102

    Fi1 = {"ao_log": [1, 3.00160],
       "pow": [0, 1], "ao_pow": [9.91243972, -6.33270087],
       "ao_exp": [0.008449, 4.6942, 3.4865, 1.6572, 1.4115],
       "titao": [648/Tc, 1957/Tc, 3895/Tc, 5705/Tc, 15080/Tc]}

    Fi2 = {"R": 8.314510,
       "ao_log": [1, 3.00088],
       "pow": [0, 1], "ao_pow": [19.597508817, -83.959667892],
       "ao_exp": [], "titao": [],
       "ao_sinh": [0.76315, 8.74432], "sinh": [4.306474465, 5.577233895],
       "ao_cosh": [0.0046, -4.46921], "cosh": [0.936220902, 5.722644361]}

    Fi3 = {"ao_log": [1, 2.5998324],
       "pow": [0, -1./3, -2./3, -1],
       "ao_pow": [-10.413865, -3.3854083, 1.6900979, -0.3911541],
       "ao_exp": [4.7206715], "titao": [10.543907]}

    CP4 = {"ao": 0.15438149595e2,
       "an": [-0.18044750507e7, 0.77426666393e5, -0.13241658754e4,
              -0.51479005257e-1, 0.10809172196e-3, -0.65501783437e-7],
       "pow": [-3, -2, -1, 1, 2, 3],
       "ao_exp": [-0.67490056171e1], "exp": [3000]}




    GERG = {
    "__type__": "Helmholtz",
    "__name__": "Helmholtz equation of state for methane of Kunz and "
                "Wagner (2004).",
    "__doi__": {"autor": "Kunz, O., Wagner, W.",
                "title": "The GERG-2008 Wide-Range Equation of State for "
                         "Natural Gases and Other Mixtures: An Expansion "
                         "of GERG-2004",
                "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                "doi":  "10.1021/je300655b"},
    "R": 8.314472,
    "cp": Fi2,
    "ref": "OTO",

    "Tmin": 90.6941, "Tmax": 625.0, "Pmax": 1000000.0, "rhomax": 40.072,

    "nr1":  [0.57335704239162, -0.16760687523730e1, 0.23405291834916,
             -0.21947376343441, 0.16369201404128e-1, 0.15004406389280e-1],
    "d1": [1, 1, 2, 2, 4, 4],
    "t1": [0.125, 1.125, 0.375, 1.125, 0.625, 1.5],

    "nr2": [0.98990489492918e-1, 0.58382770929055, -0.7478686756039,
            0.30033302857974, 0.20985543806568, -0.18590151133061e-1,
            -0.15782558339049, 0.12716735220791, -0.32019743894346e-1,
            -0.68049729364536e-1, 0.24291412853736e-1, 0.51440451639444e-2,
            -0.019084949733532, 0.55229677241291e-2, -0.44197392976085e-2,
            0.040061416708429, -0.33752085907575e-1, -0.25127658213357e-2],
    "d2": [1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7],
    "t2": [0.625, 2.625, 2.75, 2.125, 2, 1.75, 4.5, 4.75, 5, 4, 4.5, 7.5,
           14, 11.5, 26, 28, 30, 16],
    "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6],
    "gamma2": [1]*18}