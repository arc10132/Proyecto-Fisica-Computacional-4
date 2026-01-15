# Proyecto-Fisica-Computacional-4
Este proyecto trata de Migración Orbital de Partı́culas en Discos Protoplanetarios Huecos 


Migración Orbital de Partı́culas en Discos Protoplanetarios Huecos
En este proyecto se estudiará la migración de partı́culas sólidas inmersas en un disco protoplanetario que rodea a una
estrella central. El disco presenta un hueco interno, es decir, no hay gas dentro de un radio rmin alrededor de la estrella.
El gas ejerce dos efectos principales sobre las partı́culas:
1. Su gravedad, que contribuye al potencial total del sistema.
2. Su arrastre, que tiende a acelerar o frenar a las partı́culas hasta que se adapten a la velocidad local del gas.
El objetivo es simular de forma numérica cómo estas partı́culas evolucionan radialmente en el disco debido a estas interac-
ciones.
Etapas del proyecto:
1. Cálculo del potencial gravitatorio del disco y estrella
Considerar una densidad de masa radial del disco, expresada en función del radio r =
(
0
r < rmin ,
ρdisk (r) =

r0 p
ρ0 r
r ≥ rmin ,
p
x2 + y 2 :
donde rmin define el hueco interno y p es un ı́ndice de densidad radial.
La estrella central se modela como un núcleo compacto uniforme de radio R⋆ :

 M⋆ , px2 + y 2 ≤ R
⋆
ρ⋆ (x, y) = πR⋆2 p

2
2
x + y > R⋆
0,
Resolver la ecuación de Poisson en 2D cartesianas considerando ambas contribuciones:
∇2 Φ(x, y) = 4πG [ρ⋆ (x, y) + ρdisk (x, y)]
usando diferencias finitas.
Se obtiene el potencial Φ(x, y) en toda la malla.
2. Cálculo del campo gravitatorio
El campo se obtiene mediante:

g(x, y) = −∇Φ(x, y) = −
∂Φ ∂Φ
,
∂x ∂y

Discretizar el gradiente usando diferencias finitas centradas.
Este campo incluye tanto la gravedad de la estrella central como la del disco.
3. Integración de las trayectorias de las partı́culas
Cada partı́cula sigue la segunda ley de Newton con arrastre:
dr
dv
= g(x, y) − γ(v − vgas (r)),
= v.
dt
dt
El coeficiente de arrastre γ puede ser constante.
La velocidad del gas se define como un campo azimutal en función de r:
 r α
 y x
0
ϕ̂, ϕ̂ = − ,
vgas (r) = v0
r
r r
Integrar usando métodos numéricos como Runge-Kutta 4 o Verlet .
Registrar la evolución radial y las trayectorias de las partı́culas.
4. Visualización
Mapas del potencial Φ(x, y) y del campo g(x, y).
Trayectorias de partı́culas mostrando la migración radial en el disco.
Comparación de diferentes perfiles de densidad del disco, radios del hueco interno y valores de arrastre.
Nota: Este modelo simplificado no considera la presión del gas ni efectos hidrodinámicos complejos; el objetivo es centrarse
en la interacción gravitatoria y el arrastre básico para estudiar la migración de partı́culas.
