# -*- coding: utf-8 -*-

# Tarea numérica - Ecuaciones Diferenciales Ordinarias

# Nombre: Diego Alonso Sánchez Manríquez
# RUT: 19.957.060-9

#Librerías importadas
import matplotlib.pyplot as plt #usada para graficar
import numpy as np #usada para resolver vectorialmente
from scipy.integrate import solve_ivp #usada para integrar con RKF

#%%

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  @ funciones dadt() dmdt() dsdt()

  Motivación
  - Simplificar la llamada de las funciones del lado derecho de cada
    EDO asociada al modelo simple de formación de estrellas.

  Parámetros
  - a (float): fracción de masa de gas atómico
  - m (float): fracción de masa de gas molecular
  - s (float): fracción de masa de estrellas activas 
  - cte (dict): diccionario con constantes usadas

  Funcionamiento
  - Al ingresar los parámetros, se entrega la evaluación de estos
    en la respectiva función asociada a la EDO.

  Consideraciones
  - dadt() corresponde al lado derecho de la EDO asociada a da/dt
  - dmdt() corresponde al lado derecho de la EDO asociada a dm/dt
  - dsdt() corresponde al lado derecho de la EDO asociada a ds/dt
  
  Nota: se despejó la función del lado derecho asociada a ds/dt a
  partir de la relación entregada en el enunciado. Esto se hizo
  así, ya que al usar directamente s=1-a-m se obtenían valores
  complejos por errores de cómputo.
  
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def dadt(a,m,s,cte):
    
    #Condiciones de los parámetros
    assert type(a)==float
    assert type(m)==float
    assert type(s)==float  
    assert type(cte)==dict

    #Se entrega la evaluación
    return s - a*cte['k1']*m**2

def dmdt(a,m,s,cte):
    
    #Condiciones de los parámetros
    assert type(a)==float
    assert type(m)==float
    assert type(s)==float  
    assert type(cte)==dict

    #Se entrega la evaluación
    return a*cte['k1']*m**2 - cte['k2']*s*m**cte['alpha']

def dsdt(m,s,cte):

    #Condiciones de los parámetros
    assert type(m)==float
    assert type(s)==float  
    assert type(cte)==dict

    #Se entrega la evaluación
    return -s + cte['k2']*s*m**cte['alpha']

#%%

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  @ función euler_progresivo()

  Motivación
  - Utilizar el método de Euler (progresivo) para resolver el sistema
    de EDO's asociado a las funciones a(t), m(t) y s(t).

  Parámetros
  - T (int): extremo superior del intervalo a analizar
  - dt (float): paso de tiempo (medido en millones de años)
  - cte (dict): diccionario con constantes usadas

  Funcionamiento
  - Al llamar la función, se entregan cuatro listas (t,a,m,s) que
    corresponden a la solución numérica del sistema de EDO's.
  
  Consideración
  - Como la función entrega cuatro listas, se deben "recibir" con 
    una asignación múltiple de la forma t,a,m,s=euler_progresivo().
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def euler_progresivo(T,dt,cte):

    #Condiciones de los parámetros
    assert type(T)==int
    assert type(dt)==float
    assert type(cte)==dict

    #Cantidad de puntos 
    N = int(T/dt)
    
    #Condiciones iniciales
    t0 = 0; a0 = cte['a0']; m0 = cte['m0']; s0 = 1-a0-m0

    #Creación de las listas donde se guardan las soluciones
    t = [t0]; a = [a0]; m = [m0]; s = [s0] 

    #Se aplica Euler (progresivo) guardando los valores
    for i in range(0,N):
        t += [t[i] + dt]
        a += [a[i] + dt*dadt(a[i],m[i],s[i],cte)]
        m += [m[i] + dt*dmdt(a[i],m[i],s[i],cte)]
        s += [s[i] + dt*dsdt(m[i],s[i],cte)]

    #Se entregan las soluciones al sistema de EDO's
    return t,a,m,s

#%%

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  @ función runge_kutta4()

  Motivación
  - Utilizar el método Runge-Kutta de orden 4 para resolver el
    sistema de EDO's asociado a las funciones a(t), m(t) y s(t).

  Parámetros
  - T (int): extremo superior del intervalo a analizar
  - dt (float): paso de tiempo (medido en millones de años)
  - cte (dict): diccionario con constantes usadas

  Funcionamiento
  - Al llamar la función, se entregan cuatro listas (t,a,m,s) que
    corresponden a la solución numérica del sistema de EDO's.
  
  Consideración
  - Como la función entrega cuatro listas, se deben "recibir" con 
    una asignación múltiple de la forma t,a,m,s=runge_kutta4().
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def runge_kutta4(T,dt,cte):

    #Condiciones de los parámetros
    assert type(T)==int
    assert type(dt)==float
    assert type(cte)==dict    

    #Funciones auxiliares para los pasos de Runge-Kutta 4
    
    def _dadt(dt,a,m,s): #paso para lado derecho de dadt
        a1 = dadt(a,m,s,cte)
        a2 = dadt(a+a1*dt/2,m+a1*dt/2,s+a1*dt/2,cte)
        a3 = dadt(a+a2*dt/2,m+a2*dt/2,s+a2*dt/2,cte)
        a4 = dadt(a+a3*dt/2,m+a3*dt/2,s+a3*dt/2,cte)
        return (a1+2*a2+2*a3+a4)*dt/6

    def _dmdt(dt,a,m,s): #paso para lado derecho de dmdt
        m1 = dmdt(a,m,s,cte)
        m2 = dmdt(a+m1*dt/2,m+m1*dt/2,s+m1*dt/2,cte)
        m3 = dmdt(a+m2*dt/2,m+m2*dt/2,s+m2*dt/2,cte)
        m4 = dmdt(a+m3*dt/2,m+m3*dt/2,s+m3*dt/2,cte)
        return (m1+2*m2+2*m3+m4)*dt/6

    def _dsdt(dt,m,s): #paso para lado derecho de dsdt
        s1 = dsdt(m,s,cte)
        s2 = dsdt(m+s1*dt/2,s+s1*dt/2,cte)
        s3 = dsdt(m+s2*dt/2,s+s2*dt/2,cte)
        s4 = dsdt(m+s3*dt/2,s+s3*dt/2,cte)
        return (s1+2*s2+2*s3+s4)*dt/6

    #Cantidad de puntos
    N = int(T/dt)
    
    #Condiciones iniciales
    t0 = 0; a0 = cte['a0']; m0 = cte['m0']; s0 = 1-a0-m0

    #Creación de las listas donde se guardan las soluciones
    t = [t0]; a = [a0]; m = [m0]; s = [s0] 

    #Se aplica Runge-Kutta 4 guardando los valores
    for i in range(0,N):
        
        t += [t[i] + dt]
        a += [a[i] + _dadt(dt,a[i],m[i],s[i])]
        m += [m[i] + _dmdt(dt,a[i],m[i],s[i])]
        s += [s[i] + _dsdt(dt,m[i],s[i])]
        
    #Se entregan las soluciones al sistema de EDO's
    return t,a,m,s

#%%

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  @ función solucion_RKF()

  Motivación
  - Implementar el método de Runge-Kutta-Fehlberg al modelo simple
    de formación de estrellas. Para esto se usa la función solve_ivp 
    de la librería scipy para resolver el sistema de EDO's asociado a
    las funciones a(t), m(t) y s(t).

  Parámetros
  - T (int): extremo superior del intervalo a analizar
  - dt (float): paso de tiempo (medido en millones de años)
  - cte (dict): diccionario con constantes usadas

  Funcionamiento
  - Al llamar la función, se entregan cuatro listas (t,a,m,s) que
    corresponden a la solución numérica del sistema de EDO's.
  
  Consideración
  - Como la función entrega cuatro listas, se deben "recibir" con 
    una asignación múltiple de la forma t,a,m,s=solucion_RKF().

  Nota: A diferencia de los métodos implementados anteriormente,
  este resuelve el sistema de manera vectorial para asegurar el 
  correcto funcionamiento de solve_ivp, pues trabaja con vectores.
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def solucion_RKF(T,dt,cte):
    
    #Condiciones de los parámetros
    assert type(T)==int
    assert type(dt)==float
    assert type(cte)==dict    
    
    #Función vectorial de a,m,s
    def F(dt,ams,k1,k2,alpha):
    
        #Se separan los valores
        a,m,s = ams
        
        #Se entrega el arreglo
        return [s - a*k1*m**2, #función asociada a da/dt
                a*k1*m**2 - k2*s*m**alpha, #función asociada a dm/dt
                -s + k2*s*m**alpha] #función asociada a ds/dt
    
    #Cantidad de puntos
    N = int(T/dt)
    
    #Vector de estado
    ams0 = [cte["a0"],cte["m0"],1-cte["a0"]-cte["m0"]]
    
    #Intervalos equiespaciados para la solución
    t = np.linspace(0,T,N)

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""
      Parámetros de la función solve_ivp:
      - fun: funciones del sistema expresada vectorialmente
      - t_span: intervalo de interés para la integración
      - y0: vector de estado (condiciones iniciales)
      - method: método utilizado para la resolución
      - t_eval: tiempos que se desean guardar
      - args: argumentos adicionales de fun
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    #Se aplica solve_ivp guardando los valores en el vector ams
    ams = solve_ivp(fun=F,t_span=(0,T),y0=ams0,method="RK45",\
          t_eval=t,args=(cte["k1"],cte["k2"],cte["alpha"]))
    
    #Se extraen las soluciones del vector
    t,a,m,s = ams.t,ams.y[0],ams.y[1],ams.y[2] 

    #Se entregan las soluciones al sistema de EDO's
    return t,a,m,s        

#%%

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  @ función periodo()

  Motivación
  - Encontrar el periodo límite de un sistema que tiende a tener un
    equilibrio periódico, es decir, que tiende a ser periódica.
  
  Parámetros
  - t (list): valores tomados por el tiempo
  - s (list): valores de s(t) en función del tiempo
  
  Funcionamiento
  - Al llamar la función, se entrega el último periodo encontrado.
    Esto considerando la hipótesis del enunciado, es decir, que el
    sistema tiende a un equilibrio periódico).

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def periodo(t,s):
 
    #Condiciones de los parámetros
    assert type(t)==float or int
    assert type(s)==float or int

    maximosLocales=[]
    
    #Busca los máximos locales en el intervalo estudiado
    for i in range(0,len(s)-1):

        #Se encuentra máximo local
        if s[i-1]<s[i] and s[i]>s[i+1] and round(abs(max(s)-s[i]),1)==0:
            maximosLocales+=[(t[i],s[i])] # Se agrega

    #Entrega el último intervalo suponiendo la hipótesis entregada
    return maximosLocales[-1][0]-maximosLocales[-2][0]

#%%

T = 100 #Intervalo a estudiar
dt = 0.001 #Paso del tiempo

#Diccionario con las constantes para el caso de alphas [1.3-1.9]
alpha13 = {'k1':8,'k2':15,'alpha':1.3,'a0':0.4,'m0':0.2}
alpha14 = {'k1':8,'k2':15,'alpha':1.4,'a0':0.4,'m0':0.3} 
alpha15 = {'k1':8,'k2':15,'alpha':1.5,'a0':0.4,'m0':0.3}
alpha16 = {'k1':8,'k2':15,'alpha':1.6,'a0':0.4,'m0':0.3}
alpha17 = {'k1':8,'k2':15,'alpha':1.7,'a0':0.4,'m0':0.3}
alpha18 = {'k1':8,'k2':15,'alpha':1.8,'a0':0.4,'m0':0.3}
alpha19 = {'k1':8,'k2':15,'alpha':1.9,'a0':0.4,'m0':0.3}

#Nota: dt = 0.001 para evitar divergencia en RK4 al comparar.

#%%

#Datos para parte D (Caso alpha = 1.3) 
t13_EP,a13_EP,m13_EP,s13_EP = euler_progresivo(T,dt,alpha13)
t13_RK4,a13_RK4,m13_RK4,s13_RK4 = runge_kutta4(T,dt,alpha13)
t13_RKF,a13_RKF,m13_RKF,s13_RKF = solucion_RKF(T,dt,alpha13)

#%%

#Datos para parte D (Caso alpha = 1.4) 
t14_EP,a14_EP,m14_EP,s14_EP = euler_progresivo(T,dt,alpha14)
t14_RK4,a14_RK4,m14_RK4,s14_RK4 = runge_kutta4(T,dt,alpha14)
t14_RKF,a14_RKF,m14_RKF,s14_RKF = solucion_RKF(T,dt,alpha14)

#%%

#Datos para parte D (Caso alpha = 1.5) 
t15_EP,a15_EP,m15_EP,s15_EP = euler_progresivo(T,dt,alpha15)
t15_RK4,a15_RK4,m15_RK4,s15_RK4 = runge_kutta4(T,dt,alpha15)
t15_RKF,a15_RKF,m15_RKF,s15_RKF = solucion_RKF(T,dt,alpha15)

#%%

#Datos para parte D (Caso alpha = 1.6) 
t16_EP,a16_EP,m16_EP,s16_EP = euler_progresivo(T,dt,alpha16)
t16_RK4,a16_RK4,m16_RK4,s16_RK4 = runge_kutta4(T,dt,alpha16)
t16_RKF,a16_RKF,m16_RKF,s16_RKF = solucion_RKF(T,dt,alpha16)

#%%

#Datos para parte D (Caso alpha = 1.7) 
t17_EP,a17_EP,m17_EP,s17_EP = euler_progresivo(T,dt,alpha17)
t17_RK4,a17_RK4,m17_RK4,s17_RK4 = runge_kutta4(T,dt,alpha17)
t17_RKF,a17_RKF,m17_RKF,s17_RKF = solucion_RKF(T,dt,alpha17)

#%%

#Datos para parte D (Caso alpha = 1.8) 
t18_EP,a18_EP,m18_EP,s18_EP = euler_progresivo(T,dt,alpha18)
t18_RK4,a18_RK4,m18_RK4,s18_RK4 = runge_kutta4(T,dt,alpha18)
t18_RKF,a18_RKF,m18_RKF,s18_RKF = solucion_RKF(T,dt,alpha18)

#%%

#Datos para parte D (Caso alpha = 1.9) 
t19_EP,a19_EP,m19_EP,s19_EP = euler_progresivo(T,dt,alpha19)
t19_RK4,a19_RK4,m19_RK4,s19_RK4 = runge_kutta4(T,dt,alpha19)
t19_RKF,a19_RKF,m19_RKF,s19_RKF = solucion_RKF(T,dt,alpha19)

#%%

#Lista con los valores de alpha
alphas = [alpha13['alpha'],alpha14['alpha'],alpha15['alpha'],\
          alpha16['alpha'],alpha17['alpha'],alpha18['alpha'],\
          alpha19['alpha']]

#Listas con los periodos usando Euler progresivo para cada alpha
periodos_EP = [periodo(t13_EP,s13_EP),periodo(t14_EP,s14_EP),\
               periodo(t15_EP,s15_EP),periodo(t16_EP,s16_EP),\
               periodo(t17_EP,s17_EP),periodo(t18_EP,s18_EP),\
               periodo(t19_EP,s19_EP)]

#Listas con los periodos usando Runge-Kutta 4 para cada alpha
periodos_RK4 = [periodo(t13_RK4,s13_RK4),periodo(t14_RK4,s14_RK4),\
                periodo(t15_RK4,s15_RK4),periodo(t16_RK4,s16_RK4),\
                periodo(t17_RK4,s17_RK4),periodo(t18_RK4,s18_RK4),\
                periodo(t19_RK4,s19_RK4)]    
    
#Listas con los periodos usando Runge-Kutta-Fehlberg para cada alpha
periodos_RKF = [periodo(t13_RKF,s13_RKF),periodo(t14_RKF,s14_RKF),\
                periodo(t15_RKF,s15_RKF),periodo(t16_RKF,s16_RKF),\
                periodo(t17_RKF,s17_RKF),periodo(t18_RKF,s18_RKF),\
                periodo(t19_RKF,s19_RKF)]
    
#%%

#Creación del gráfico de "periodos límites en función de alpha"
    
#Creación del lienzo
fig, ax = plt.subplots(figsize=(12,6))

#Gráficos
ax.plot(alphas,periodos_EP,label="Periodos obtenidos usando Euler progresivo")
ax.plot(alphas,periodos_RK4,label="Periodos obtenidos usando Runge-Kutta 4")
ax.plot(alphas,periodos_RKF,label="Periodos obtenidos usando Runge-Kutta-Fehlberg")

#Etiquetas
ax.set_xlabel("Valor de $\\alpha$",labelpad=20)
ax.set_ylabel("Periodo límite (millones de años)",labelpad=20)

#Título
ax.set_title("Periodo límite en función de alpha",\
             fontweight="bold", loc='center',pad=20)

#Configuraciones
ax.grid(b=True, which='major', axis='both')
ax.margins(0.1)

# Leyendas
ax.legend()

#Guardado de figura
fig.savefig('Imágenes/E Periodo en función de alpha (EP vs RK4 vs RKF).pdf')

