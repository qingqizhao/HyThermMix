import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.GERG import GERG

if __name__ == '__main__':
    T = 323
    P = 1e7  # Pa
    componente = [2, 12]  # indices for CO2, H2O for example
    fraccion = [0.4, 0.6]
    aire = GERG(T=T, P=P, componente=componente, fraccion=fraccion)
    print("T:", aire.T)
    print("P:", aire.P)
    print("rho:", aire.rho)
    print("h:", aire.h)
    print("s:", aire.s)
    print("Fugacity Coefficients:", aire.FI)
