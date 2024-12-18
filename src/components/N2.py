from math import exp
# from unittest import TestCase

from scipy.constants import Boltzmann, Avogadro


# from lib import unidades
# from lib.meos import MEoS


class N2(object):
    """Multiparamente equation of state for nitrogen"""
    name = "nitrogen"
    CASNumber = "7727-37-9"
    formula = "N2"
    synonym = "R-728"
    _refPropName = "NITROGEN"
    _coolPropName = "Nitrogen"
    # rhoc = unidades.Density(313.299958972)
    A = -8.08926E-02
    B = -9.03072E-01
    C = 2.96500E-02
    D = 1.94214E-03
    rhocL = 11.1839
    rhoc = 313.299958972
    # Tc = unidades.Temperature(126.192)
    Tc = 126.192
    # Pc = unidades.Pressure(3395.8, "kPa")
    Pc = 3395.8 * 1e3  # "kPa"
    M = 28.01348  # g/mol
    # Tt = unidades.Temperature(63.151)
    Tt = 63.151
    # Tb = unidades.Temperature(77.355)
    Tb = 77.355
    f_acent = 0.0372
    # momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 46
    # _Tr = unidades.Temperature(122.520245)
    _Tr = 122.520245
    # _rhor = unidades.Density(316.134310)
    _rhor = 316.134310
    _w = 0.043553140

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [-12.76952708, -0.00784163, -1.934819e-4,
                      -1.247742e-5, 6.678326e-8],
           "ao_exp": [1.012941],
           "titao": [26.65788]}

    Fi2 = {"ao_log": [1, 2.50031],
           "pow": [0, 1],
           "ao_pow": [11.083407489, -22.202102428],
           "ao_sinh": [0.13732, 0.90066], "sinh": [5.25182262, 13.788988208],
           "ao_cosh": [-0.1466], "cosh": [-5.393067706]}

    CP1 = {"ao": 3.50404228308756,
           "an": [-0.735210401157252e3, 0.342239980411978e2, -0.55764828456762,
                  -1.73390185081005e-5, 1.74650849766463e-8,
                  -3.56892033544348e-12],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [1.00538722808834], "exp": [3353.4061]}

    CP2 = {"ao": 3.50418363823,
           "an": [-0.837079888737e3, 0.379147114487e2, -0.601737844275,
                  -0.874955653028e-5, 0.148958507239e-7, -0.256370354277e-11],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [1.00773735767], "exp": [3353.4061]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nitrogen of Kunz and "
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

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 2200000.0, "rhomax": 53.15,

        "nr1": [0.59889711801201, -0.16941557480731e1, 0.24579736191718,
                -0.23722456755175, 0.17954918715141e-1, 0.14592875720215e-1],
        "d1": [1, 1, 2, 2, 4, 4],
        "t1": [0.125, 1.125, 0.375, 1.125, 0.625, 1.5],

        "nr2": [0.10008065936206, 0.73157115385532, -0.88372272336366,
                0.31887660246708, 0.20766491728799, -0.19379315454158e-1,
                -0.16936641554983, 0.13546846041701, -0.33066712095307e-1,
                -0.60690817018557e-1, 0.12797548292871e-1, 0.58743664107299e-2,
                -0.018451951971969, 0.47226622042472e-2, -0.52024079680599e-2,
                0.043563505956635, -0.36251690750939e-1, -0.28974026866543e-2],
        "d2": [1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7],
        "t2": [0.625, 2.625, 2.75, 2.125, 2, 1.75, 4.5, 4.75, 5, 4, 4.5, 7.5,
               14, 11.5, 26, 28, 30, 16],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6],
        "gamma2": [1] * 18,

        "nr3": [],
        "nr4": []}

    eq = GERG