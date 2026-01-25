"""
Created on 2026-01-23 

"""
# Importamos las librerias necesarias para los calculos dentro del modulo

import numpy as np
import scipy.constants as consts


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

def star_density(x , y , Mass, Radius):
    """"
    Calcula el perfil de densidad de una estrella.
    Args:
        x (Array):  Eje x de la malla
        y (Array):  Eje y de la malla      
        Mass (Float): Masa de la estrella, idealmente con unidades de kg      
        Radius (Floar): Radio de la estrella, idealmente con unidades de au     

    Returns:
        Density(Array(): Densidad superficial de la estrella en cada punto de la malla, condicionada al tamaño del radio 
    """
    x=x
    y=y
    r= np.sqrt(x**2 + y**2)
    star_density = np.where(r <= Radius, Mass / (np.pi * Radius**2) , 0 )

    return star_density


def Poisson_equation(x , y ,dx, star_density, density_profile):
    """
    """
    phi = np.zeros_like(x)

    for i in range(1, x.shape[0]-1):
        for j in range(1, y.shape[1]-1):
            phi[i,j]= ((phi[i,j+1]+phi[i+1,j]+phi[i,j-1]+phi[i-1,j])/4)- np.pi* consts.G * (star_density[i,j] +density_profile[i,j] *dx**2)
            
    return phi 
          