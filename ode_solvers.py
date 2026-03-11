# -*- coding: utf-8 -*-
"""


@author: [Omar Hernandez]
@mail: []

Description: Módulo de solucionadores de EDOs
"""

# importamos librerías
import numpy as np
import matplotlib.pyplot as plt

from typing import Callable, Union, Tuple

# definición de métodos numéricos para resolver EDOs
def euler_method_second_order_2D(
    f : Callable[[float, np.ndarray, np.ndarray], np.ndarray],
    t0: float,
    r0: Union[Tuple[float, float], np.ndarray],
    v0: Union[Tuple[float, float], np.ndarray],
    tf: float,
    h : float = 1e-5 
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Resuelve una ecuación diferencial de segundo orden en 2D de la forma:
        v' = f(t, y, v)
        y' = v
    utilizando el método de Euler.

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

    # Convertir r0 y v0 a np.ndarray si son tuplas
    r0 = np.array(r0, dtype=float)
    v0 = np.array(v0, dtype=float)
    
    # número de pasos
    n_steps = int((tf - t0) / h) + 1
    t = np.linspace(t0, tf, n_steps)
    
    # inicializar arrays para posición y velocidad
    r = np.zeros((n_steps, 2))
    v = np.zeros((n_steps, 2))
    
    # condiciones iniciales
    r[0] = r0
    v[0] = v0
    
    # método de Euler este es el ciclo que cambia por cada situacion o ejercicio
    for n in range(n_steps - 1):
        i = n
    
        a = f(t[i], r[i], v[i])  # aceleración
        r[i+1] = r[i] + h * v[i]
        v[i+1] = v[i] + h * a
    
    return t, r, v
    

def rk2_method_second_order_2D(
    
    f : Callable[[float, np.ndarray, np.ndarray], np.ndarray],
    t0: float,
    r0: Union[Tuple[float, float], np.ndarray],
    v0: Union[Tuple[float, float], np.ndarray],
    tf: float,
    h : float = 1e-5 
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Resuelve una ecuación diferencial de segundo orden en 2D de la forma:
        v' = f(t, y, v)
        y' = v
    utilizando el método de Euler.

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

    # Convertir r0 y v0 a np.ndarray si son tuplas
    r0 = np.array(r0, dtype=float)
    v0 = np.array(v0, dtype=float)
    
    # número de pasos
    n_steps = int((tf - t0) / h) + 1
    t = np.linspace(t0, tf, n_steps)
    
    # inicializar arrays para posición y velocidad
    r = np.zeros((n_steps, 2))
    v = np.zeros((n_steps, 2))
    
    # condiciones iniciales
    r[0] = r0
    v[0] = v0
    
    # método de Euler este es el ciclo que cambia por cada situacion o ejercicio
    for n in range(n_steps - 1):
        k1 = h*f(t[n],r[n],v[n])
        m1= h*v[n]
        
        k2 = h* f(t[n]+h/2,r[n]+m1/2,v[n]+k1/2)
        m2 = h*(v[n]+k1/2)
        
        v[n+1] = v[n] + k2
        r[n+1] = r[n] + m2
    return t, r, v

def rk4_method_second_order_2D(f : Callable[[float, np.ndarray, np.ndarray], np.ndarray],
    t0: float,
    r0: Union[Tuple[float, float], np.ndarray],
    v0: Union[Tuple[float, float], np.ndarray],
    tf: float,
    h : float = 1e-5 
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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

    # Convertir r0 y v0 a np.ndarray si son tuplas
    r0 = np.array(r0, dtype=float)
    v0 = np.array(v0, dtype=float)
    
    # número de pasos
    n_steps = int((tf - t0) / h) + 1
    t = np.linspace(t0, tf, n_steps)
    
    # inicializar arrays para posición y velocidad
    r = np.zeros((n_steps, 2))
    v = np.zeros((n_steps, 2))
    
    # condiciones iniciales
    r[0] = r0
    v[0] = v0
    
    # método de RK4 este es el ciclo que cambia por cada situacion o ejercicio
    for n in range(n_steps - 1):
        k1 = h*f(t[n],r[n],v[n])
        m1= h*v[n]
        
        k2 = h* f(t[n]+h/2,r[n]+m1/2,v[n]+k1/2)
        m2 = h*(v[n]+k1/2)
        
        k3 = h* f(t[n]+h/2,r[n]+m2/2,v[n]+k2/2)
        m3 = h*(v[n]+k2/2)
        
        k4= h* f(t[n]+h,r[n]+m3,v[n]+k3)
        m4=h*(v[n]+k3)
        
        
        v[n+1] = v[n] + (1/6)*(k1+2*k2+2*k3+k4)
        r[n+1] = r[n] + (1/6)*(m1+2*m2+2*m3+m4)
        
        
    return t, r, v

def verlet_method():
    pass