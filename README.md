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
\dfrac{M_\star}{\pi R_\star^2} & \sqrt{x^2 + y^2} \leq R_\star,
0 & \sqrt{x^2 + y^2} > R_\star,
\end{cases}
$$

donde:
- $M_\star$ es la masa total de la estrella
- $R_\star$ es el radio estelar
- La densidad es constante dentro del radio $R_\star$
- Fuera de $R_\star$, la densidad es nula



