# Migración Orbital de Partículas en Discos Protoplanetarios 

En este proyecto se estudiará la migración de partículas sólidas inmersas en un disco protoplanetario que rodea a una estrella central. El disco presenta un **hueco interno**, es decir, no hay gas dentro de un radio r_min alrededor de la estrella. El gas ejerce dos efectos principales sobre partículas:

1.- Su **gravedad**, que contribuye al potencial total del sistema.

2.- Su **arrastre**, que tiende a acelerar o frenar a las partículas evolucionan radialmente en el disco debido a estas interacciones.

## Partes del proyecto:

Se calcular el potencial grtivatoria del disco y de la estrella, con las siguientes funciones:

La densidad del disco protoplanetario se calcula como:

$$
\rho_{\text{disk}}(r) = 
\begin{cases} 
0 & r < r_{\text{min}}, \\
\rho_0 \left( \dfrac{r_0}{r} \right)^p & r \geq r_{\text{min}},
\end{cases}
$$

donde:
- $\rho_0$ es la densidad de referencia
- $r_0$ es el radio de referencia  
- $p$ es el índice de densidad
- $r_{\text{min}}$ es el radio del hueco interno


La estrella central se modela como un núcleo compacto uniforme de radio $R_\star$:

$$
\rho_\star(x, y) =
\begin{cases} 
\dfrac{M_\star}{\pi R_\star^2} & \sqrt{x^2 + y^2} \leq R_\star, \\
0 & \sqrt{x^2 + y^2} > R_\star,
\end{cases}
$$

donde:
- $M_\star$ es la masa total de la estrella
- $R_\star$ es el radio estelar
- La densidad es constante dentro del radio $R_\star$
- Fuera de $R_\star$, la densidad es nula




Ecuacion de Poisson:
El potencial gravitacional $\Phi(x,y)$ satisface la ecuación de Poisson:

$$
\nabla^2 \Phi(x,y) = 4\pi G \left[\rho_\star(x,y) + \rho_{\text{disk}}(x,y)\right]
$$

donde:
- $\nabla^2$ es el operador Laplaciano en 2D
- $G$ es la constante gravitacional
- $\rho_\star$ es la densidad de la estrella
- $\rho_{\text{disk}}$ es la densidad del disco


El campo gravitacional se obtiene como el gradiente negativo del potencial:

$$
\mathbf{g}(x,y) = -\nabla \Phi(x,y) = -\left(\frac{\partial\Phi}{\partial x}, \frac{\partial\Phi}{\partial y}\right)
$$

donde:
- $\nabla$ es el operador gradiente
- $\mathbf{g}$ es la aceleración gravitacional


La dinámica de las partículas está gobernada por:

$$
\begin{aligned}
\frac{d\mathbf{v}}{dt} &= \mathbf{g}(x, y) - \gamma (\mathbf{v} - \mathbf{v}_{\text{gas}}(r)) \\
\frac{d\mathbf{r}}{dt} &= \mathbf{v}
\end{aligned}
$$

donde:
- $\mathbf{v}$ es la velocidad de la partícula
- $\gamma$ es el coeficiente de arrastre (fricción)
- $\mathbf{v}_{\text{gas}}$ es la velocidad del gas


El gas del disco rota con un perfil de velocidad tipo Kepler:

$$
\mathbf{v}_{\text{gas}}(r) = v_0 \left( \frac{r_0}{r} \right)^\alpha \hat{\phi}, \quad 
\hat{\phi} = \left( -\frac{y}{r}, \frac{x}{r} \right)
$$

donde:
- $v_0$ es la velocidad en el radio de referencia $r_0$
- $\alpha$ es el índice de velocidad (típicamente $\alpha = 0.5$ para Kepleriano)
- $\hat{\phi}$ es el vector unitario azimutal
- $r = \sqrt{x^2 + y^2}$ es la distancia radial





