# Migración Orbital de Partículas en Discos Protoplanetarios

Este proyecto se realizó como trabajo final para el ramo **Física Computacional IV** de la carrera **Astrofísica con mención en Ciencia de Datos** de la Universidad de Santiago de Chile.  
El objetivo es estudiar la migración de **partículas sólidas** inmersas en un disco protoplanetario que rodea a una estrella central, considerando efectos gravitacionales y de arrastre del gas.

El disco presenta un **hueco interno**, es decir, no hay gas dentro de un radio `r_min` alrededor de la estrella. El gas ejerce dos efectos principales sobre las partículas:

1. **Gravedad**, pertenecientes a la estrella y al disco, que contribuye al potencial total del sistema.  
2. **Arrastre** del gas, que tiende a acelerar o frenar a las partículas, provocando migración radial en el disco.

---

## Integrantes

- Harry González Garrido  
- Tomás Inzulza Mejías  

---

## Modelo Físico

### Estrella Central

Se modela como un **núcleo compacto uniforme** de radio $R_\star$:

$$
\rho_\star(x, y) =
\begin{cases} 
\dfrac{M_\star}{\pi R_\star^2} & \sqrt{x^2 + y^2} \leq R_\star, \\
0 & \sqrt{x^2 + y^2} > R_\star,
\end{cases}
$$

- $M_\star$ es la masa total de la estrella  
- $R_\star$ es el radio estelar  

### Disco Protoplanetario

La densidad del disco se calcula como:

$$
\rho_{\text{disk}}(r) = 
\begin{cases} 
0 & r < r_{\text{min}}, \\
\rho_0 \left( \dfrac{r_0}{r} \right)^p & r \geq r_{\text{min}},
\end{cases}
$$

- $\rho_0$: densidad de referencia  
- $r_0$: radio de referencia  
- $p$: índice de densidad  
- $r_{\text{min}}$: radio del hueco interno  

### Campo Gravitacional

Se resuelve la **ecuación de Poisson** para obtener el potencial gravitatorio:

$$
\nabla^2 \Phi(x,y) = 4\pi G \left[\rho_\star(x,y) + \rho_{\text{disk}}(x,y)\right]
$$

Y el campo gravitacional se obtiene como el gradiente negativo del potencial:

$$
\mathbf{g}(x,y) = -\nabla \Phi(x,y) = -\left(\frac{\partial\Phi}{\partial x}, \frac{\partial\Phi}{\partial y}\right)
$$

### Gas del Disco

El gas rota con un perfil tipo Kepler:

$$
\mathbf{v}_{\text{gas}}(r) = v_0 \left( \frac{r_0}{r} \right)^\alpha \hat{\phi}, \quad 
\hat{\phi} = \left( -\frac{y}{r}, \frac{x}{r} \right)
$$

- $v_0$: velocidad en el radio de referencia $r_0$  
- $\alpha$: índice de velocidad  
- $\hat{\phi}$: vector unitario azimutal  
- $r = \sqrt{x^2 + y^2}$: distancia radial  

### Ecuaciones de Movimiento

La dinámica de las partículas se describe con la **segunda ley de Newton**, considerando arrastre:

$$
\begin{aligned}
\frac{d\mathbf{v}}{dt} &= \mathbf{g}(x, y) - \gamma (\mathbf{v} - \mathbf{v}_{\text{gas}}(r)) \\
\frac{d\mathbf{r}}{dt} &= \mathbf{v}
\end{aligned}
$$

- $\mathbf{v}$: velocidad de la partícula  
- $\gamma$: coeficiente de arrastre  
- $\mathbf{v}_{\text{gas}}$: velocidad del gas  

La integración se realiza con **Runge-Kutta de cuarto orden (RK4)**.

---

## Solución Numérica

- Discretización del potencial gravitacional mediante la ecuación de Poisson en una malla 2D  
- Interpolación del campo gravitacional con **RegularGridInterpolator** de SciPy  
- Integración de trayectorias de partículas usando RK4  
- Visualización de densidades, potencial, campo gravitacional y trayectorias  

---

## Resultados

- Se obtienen mapas de densidad para la estrella y el disco.  
- El **potencial gravitatorio** es mínimo en la región central, generando un gradiente radial hacia la estrella.  
- El **campo gravitacional** dirige a las partículas hacia el centro, combinando con el arrastre del gas para producir migración radial.  
- Se analizaron múltiples casos de trayectorias:  
  - Partículas lanzadas desde distintos radios  
  - Con velocidad inicial horizontal  
  - Con intento de escape vertical o diagonal  
  - Lanzadas directamente hacia el sistema, mostrando órbitas caóticas  

---

## Conclusión

- El modelo es una **aproximación simplificada** de un sistema estrella-disco con hueco interno.  
- Permite estudiar la **dinámica de partículas sólidas**, considerando gravedad y arrastre.  
- Aunque no incluye interacciones N-cuerpos ni efectos hidrodinámicos, es útil para explorar procesos de migración orbital y evolución de partículas en discos protoplanetarios.  

---

## Características del Proyecto

- Modelado de densidad de disco y estrella  
- Resolución del potencial gravitacional mediante Poisson  
- Campo gravitacional interpolado con SciPy  
- Integración de trayectorias con RK4 incluyendo arrastre  
- Visualización de mapas y trayectorias  
