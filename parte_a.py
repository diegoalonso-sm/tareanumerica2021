# -*- coding: utf-8 -*-

# Tarea numérica - Ecuaciones Diferenciales Ordinarias

# Nombre: Diego Alonso Sánchez Manríquez
# RUT: 19.957.060-9

#Librerías importadas
import matplotlib.pyplot as plt #usada para graficar

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
  @ función graficarA()

  Motivación
  - Simplificar la creación de los gráficos de la parte A. 

  Parámetros
  - t (list): valores tomados por el tiempo
  - a (list): valores de a(t) en función del tiempo
  - m (list): valores de m(t) en función del tiempo
  - s (list): valores de s(t) en función del tiempo
  - caso (str): nombre para la variante del gráfico 
  
  Funcionamiento
  - Al llamar la función, se utilizan las listas t,a,m,s para 
    graficar las curvas (t,a(t)),(t,m(t)),(t,s(t)) y guardar el
    gráfico generado en la carpeta "Imágenes" en el directorio 
    del código.
  
  Consideración
  - Internamente, la función toma el string de la variable "caso" 
    para el título y nombre con el que se guarda en "Imágenes".

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def graficarA(t,a,m,s,caso):
    
    #Condiciones de los parámetros
    assert type(t)==list
    assert type(a)==list
    assert type(m)==list
    assert type(s)==list
    assert type(caso)==str
    
    #Creación del lienzo
    fig, ax = plt.subplots(figsize=(12,6))
    
    #Gráficos
    ax.plot(t,a,label="Fracción de masa de gas atómico $a(t)$")
    ax.plot(t,m,label="Fracción de masa de gas molecular $m(s)$")
    ax.plot(t,s,label="Fracción de masa de estrellas activas $s(t)$")
    
    #Etiquetas
    ax.set_xlabel("Tiempo (millones de años)",labelpad=20)
    ax.set_ylabel("Fracción de masa",labelpad=20)
    
    #Título
    ax.set_title(f"Fracciones de masa del sistema a lo largo del tiempo ({caso})",\
                  fontweight="bold", loc='center', pad=20)
    
    #Leyendas
    ax.legend(loc="center right")
    
    #Configuraciones
    ax.grid(b=True, which='major', axis='both')
    ax.margins(0.1)
    
    #Guardado de figura
    plt.savefig(f'Imágenes/A ({caso}).pdf')

#%%

T = 100 #Intervalo a estudiar
dt = 0.1 #Paso del tiempo

#Diccionario con las constantes usadas para cada caso
caso1 = {'k1':10,'k2':10,'alpha':1.0,'a0':0.15,'m0':0.15}
caso2 = {'k1': 8,'k2':15,'alpha':1.2,'a0':0.40,'m0':0.30} 
caso3 = {'k1': 8,'k2':15,'alpha':1.5,'a0':0.40,'m0':0.30}
caso4 = {'k1': 8,'k2':15,'alpha':1.9,'a0':0.40,'m0':0.30}
caso5 = {'k1': 8,'k2':15,'alpha':2.0,'a0':0.40,'m0':0.30}
caso6 = {'k1': 8,'k2':15,'alpha':2.1,'a0':0.40,'m0':0.30}

#Nota: dt=0.1 pues el método no diverge.

#%%

#Datos para parte A (Caso 1) 
t1,a1,m1,s1 = euler_progresivo(T,dt,caso1)

#Gráfico
graficarA(t1,a1,m1,s1,"Caso 1")

#%%

#Datos para parte A (Caso 2) 
t2,a2,m2,s2 = euler_progresivo(T,dt,caso2)

#Gráfico
graficarA(t2,a2,m2,s2,"Caso 2")

#%%

#Datos para parte A (Caso 3) 
t3,a3,m3,s3 = euler_progresivo(T,dt,caso3)

#Gráfico
graficarA(t3,a3,m3,s3,"Caso 3")

#%%

#Datos para parte A (Caso 4) 
t4,a4,m4,s4 = euler_progresivo(T,dt,caso4)

#Gráfico
graficarA(t4,a4,m4,s4,"Caso 4")

#%%

#Datos para parte A (Caso 5) 
t5,a5,m5,s5 = euler_progresivo(T,dt,caso5)

#Gráfico
graficarA(t5,a5,m5,s5,"Caso 5")

#%%

#Datos para parte A (Caso 6) 
t6,a6,m6,s6 = euler_progresivo(T,dt,caso6)

#Gráfico
graficarA(t6,a6,m6,s6,"Caso 6")

