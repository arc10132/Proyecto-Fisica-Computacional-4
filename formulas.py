"""
Created on 2026-01-23

Authors:
    Harry González Garrido
    Tomás Inzulza Mejías

Emails:
    harry.gonzalez@usach.cl
    tomas.inzulza@usach.cl

Description:
    Módulo de fórmulas para el proyecto del ramo "Física Computacional IV".
    Incluye funciones para:

    - Perfil de densidad del disco protoplanetario
    - Densidad de una estrella
    - Resolver la ecuación de Poisson para el potencial gravitacional
    - Calcular el campo gravitacional a partir del potencial
    - Método RK4 para ecuaciones diferenciales de segundo orden en 2D
"""

# Importamos las librerias necesarias para los calculos dentro del modulo
import numpy as np
import matplotlib.pyplot as plt
from typing import Callable


# Definimos la funcion para calcular el perfil de densidad del disco protoplanetario
def density_profile(r, rho_0, r_0, p, r_min):
    """
    Calcula el perfil de densidad del disco protoplanetario.
    Args:
        r (array): Radio desde el centro del disco.
        rho_0 (float): Densidad del radio r_0.
        r_0 (float): Radio de referencia.
        p (float): Indice de densidad.
        r_min (float): Radio del hueco interno.

    Returns:
        array: Densidad del disco.
    """
    density = np.where(r >= r_min, rho_0 * (r / r_0) ** (-p), 0)

    return density

# Definimos la funcion para calcular el perfil de densidad de una estrella
def star_density(x , y , Mass, Radius):
    """"
    Calcula el perfil de densidad de una estrella.
    Args:
        x (Array):  Eje x de la malla
        y (Array):  Eje y de la malla      
        Mass (Float): Masa de la estrella, idealmente con unidades de kg      
        Radius (Floar): Radio de la estrella, idealmente con unidades de au     

    Returns:
        array: Densidad superficial de la estrella en cada punto de la malla, condicionada al tamaño del radio 
    """
    x=x
    y=y
    r= np.sqrt(x**2 + y**2)
    star_density = np.where(r <= Radius, Mass / (np.pi * Radius**2) , 0 )

    return star_density

# Definimos la funcion para resolver la ecuacion de Poisson y obtener el potencial gravitacional
def poisson_equation(dx, star_density, disk_density, tol=1e-6, max_iter=50000):
    """
    Resuelve la ecuación de Poisson para obtener el potencial gravitacional.
    Args:
        dx (float): Paso espacial.
        star_density (array): Densidad de la estrella.
        disk_density (array): Densidad del disco protoplanetario.
        tol (float): Tolerancia para la convergencia.
        max_iter (int): Número máximo de iteraciones.

    Returns:
        array: Potencial gravitacional.
    """
    G = 1 # Para simplificar la integración numérica
    rho = disk_density + star_density

    phi = np.zeros_like(rho)

    for k in range(max_iter):
        phi_new = np.zeros_like(phi)
        
        phi_new[1:-1, 1:-1] = (
            (phi[2:, 1:-1] + phi[:-2, 1:-1] +
             phi[1:-1, 2:] + phi[1:-1, :-2]) / 4
            - (np.pi * G * rho[1:-1, 1:-1] * dx**2)
        )

        num = np.max(np.abs(phi_new - phi))
        den = np.max(np.abs(phi_new)) + 1e-14
        epsilon = num / den

        if epsilon < tol:
            print(f"Jacobi convergió en {k} iteraciones")
            return phi_new

        phi = phi_new.copy()

    return phi

          
# Definimos la funcion para calcular el campo gravitacional a partir del gradiente del potencial
def campo_gravitacional(Gradient_potential, dx, dy):
    """
    Calcula el campo gravitacional a partir del gradiente del potencial.
    Args:
        Gradient_potential (array): Gradiente del potencial gravitacional.
        dx (float): Paso espacial en la dirección x.
        dy (float): Paso espacial en la dirección y.
    Returns:
        gx (array): Componente x del campo gravitacional.
        gy (array): Componente y del campo gravitacional.
    """
    Nx, Ny = Gradient_potential.shape
    gx = np.zeros(shape=(Nx-2,Ny-2))
    gy = np.zeros_like(gx)


    for i in range(1, Nx-2):
        for j in range(1, Ny-2):
            # gx es la derivada en X (cambian las filas 'i')
            gx[i,j] = -(Gradient_potential[i+1,j] - Gradient_potential[i-1,j])/(2*dx) 
            # gy es la derivada en Y (cambian las columnas 'j')
            gy[i,j] = -(Gradient_potential[i,j+1] - Gradient_potential[i,j-1])/(2*dy) 

    return gx, gy


# Definimos la funcion para resolver una ecuacion diferencial de segundo orden en 2D utilizando el metodo de RK4
def rk4_method_second_order_2D(f: Callable, t0: float, r0: tuple, v0: tuple, tf: float, h: float = 0.0001) -> tuple:
    """
    Resuelve una ecuación diferencial de segundo orden en 2D de la forma:
        v' = f(t, y, v)
        y' = v
    utilizando el método de rk4.

    Args:
        f (Callable): función que define la ecuación diferencial de segundo orden.
                      Debe retornar la aceleración como un np.ndarray (2,).
        t0 (float): tiempo inicial.
        r0 (tuple o np.ndarray): posición inicial (x0, y0).
        v0 (tuple o np.ndarray): velocidad inicial (vx0, vy0).
        tf (float): tiempo final.
        h (float, opcional): tamaño de paso de tiempo. Por defecto: 1e-5.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: vectores de tiempo, posiciones y velocidades.
            - t: np.ndarray de forma (N,)
            - y: np.ndarray de forma (N, 2)
            - v: np.ndarray de forma (N, 2)
    """
    # Inicializamos vectores para guardar la solución
    t = np.arange(start=t0, step=h, stop=tf + h)
    r = np.zeros(shape=(len(t), 2))
    v = np.zeros_like(r)

    # Fijamos las condiciones iniciales
    r[0] = r0
    v[0] = v0

    # Iteramos
    for n in range(len(t) - 1):

        # Cálculo de k1 y m1
        k1 = h * f(t[n], r[n], v[n])
        m1 = h * v[n]

        # Cálculo de k2 y m2
        k2 = h * f(t[n] + h/2, r[n] + m1/2, v[n] + k1/2)
        m2 = h * (v[n] + k1/2)

        # Cálculo de k3 y m3
        k3 = h * f(t[n] + h/2, r[n] + m2/2, v[n] + k2/2)
        m3 = h * (v[n] + k2/2)

        # Cálculo de k4 y m4
        k4 = h * f(t[n] + h, r[n] + m3, v[n] + k3)
        m4 = h * (v[n] + k3)

        # Actualizamos los valores de velocidad y posición
        v[n+1] = v[n] + (k1 + 2*k2 + 2*k3 + k4) / 6
        r[n+1] = r[n] + (m1 + 2*m2 + 2*m3 + m4) / 6

    return t, r, v
