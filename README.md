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

