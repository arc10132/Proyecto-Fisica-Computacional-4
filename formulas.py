"""
Created on 2026-01-23 

"""
# Importamos las librerias necesarias para los calculos dentro del modulo

import numpy as np

# Definimos la funcion para calcular el perfil de densidad del disco protoplanetario
def density_profile(r, rho_0, r_0, p, r_min):
    """
    Calcula el perfil de densidad del disco protoplanetario.
    Args:
        r (float): Radio desde el centro del disco.
        rho_0 (float): Densidad del radio r_0.
        r_0 (float): Radio de referencia.
        p (float): Indice de densidad.
        r_min (float): Radio del hueco interno.

    Returns:
        array: Densidad del disco.
    """
    density = np.where(r >= r_min, rho_0 * (r / r_0) ** (-p), 0)

    return density
