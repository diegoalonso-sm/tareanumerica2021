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

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  @ función graficarC()

  Motivación
  - Simplificar la creación de los gráficos de la parte C.
  - Obtener los gráficos para analizar la tendencia de sus "periodos".

  Parámetros
  - t (list): valores tomados por el tiempo
  - s (list): valores de s(t) en función del tiempo
  - caso (str): nombre para la variante del gráfico 

  Funcionamiento
  - Al llamar la función, se utilizan las listas t,s para graficar 
    (t,s(t)) y guardar el gráfico generado en la carpeta "Imágenes"
    en el directorio del código.
  
  Consideración
  - Internamente, la función toma el string de la variable "caso" 
    para el título y nombre con el que se guarda en "Imágenes".

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def graficarC(t,s,caso):
    
    #Condiciones de los parámetros
    assert type(t)==list
    assert type(s)==list
    assert type(caso)==str
    
    #Creación del lienzo
    fig, ax = plt.subplots(figsize=(12,6))
    
    #Gráficos
    ax.plot(t,s)
    
    #Etiquetas
    ax.set_xlabel("Tiempo (millones de años)",labelpad=20)
    ax.set_ylabel("Fracción de masa de estrellas activas s(t)",labelpad=20)
    
    #Título
    ax.set_title(f"Fracción de masa de estrellas activas s(t) a lo largo del tiempo (alpha = {caso})",\
                 fontweight="bold", loc='center',pad=20)
    
    #Configuraciones
    ax.grid(b=True, which='major', axis='both')
    ax.margins(0.1)
    
    #Guardado de figura
    fig.savefig(f'Imágenes/C (alpha = {caso}).pdf')


#%%

T = 200 #Intervalo a estudiar
dt = 0.001 #Paso del tiempo

#Diccionario con las constantes para el caso de alphas [1.3-1.9]
alpha13 = {'k1':8,'k2':15,'alpha':1.30,'a0':0.4,'m0':0.2}
alpha14 = {'k1':8,'k2':15,'alpha':1.40,'a0':0.4,'m0':0.3} 
alpha15 = {'k1':8,'k2':15,'alpha':1.50,'a0':0.4,'m0':0.3}
alpha16 = {'k1':8,'k2':15,'alpha':1.60,'a0':0.4,'m0':0.3}
alpha17 = {'k1':8,'k2':15,'alpha':1.70,'a0':0.4,'m0':0.3}
alpha18 = {'k1':8,'k2':15,'alpha':1.80,'a0':0.4,'m0':0.3}
alpha19 = {'k1':8,'k2':15,'alpha':1.90,'a0':0.4,'m0':0.3}

#Nota: se expande el intervalo a 200 estudiar para explorar sus valores.
#Nota: dt = 0.001 para evitar divergencia en RK4 al comparar.

#%%
 
#Datos para parte C (Caso alpha = 1.3) 
t13_EP,a13_EP,m13_EP,s13_EP = euler_progresivo(T,dt,alpha13)

#Gráfico
graficarC(t13_EP,s13_EP,"1.3")

#%%

#Datos para parte C (Caso alpha = 1.4) 
t14_EP,a14_EP,m14_EP,s14_EP = euler_progresivo(T,dt,alpha14)

#Gráfico
graficarC(t14_EP,s14_EP,"1.4")

#%%

#Datos para parte C (Caso alpha = 1.5) 
t15_EP,a15_EP,m15_EP,s15_EP = euler_progresivo(T,dt,alpha15)

#Gráfico
graficarC(t15_EP,s15_EP,"1.5")

#%%

#Datos para parte C (Caso alpha = 1.6) 
t16_EP,a16_EP,m16_EP,s16_EP = euler_progresivo(T,dt,alpha16)

#Gráfico
graficarC(t16_EP,s16_EP,"1.6")

#%%

#Datos para parte C (Caso alpha = 1.7) 
t17_EP,a17_EP,m17_EP,s17_EP = euler_progresivo(T,dt,alpha17)

#Gráfico
graficarC(t17_EP,s17_EP,"1.7")

#%%

#Datos para parte C (Caso alpha = 1.8) 
t18_EP,a18_EP,m18_EP,s18_EP = euler_progresivo(T,dt,alpha18)

#Gráfico
graficarC(t18_EP,s18_EP,"1.8")

#%%

#Datos para parte C (Caso alpha = 1.9) 
t19_EP,a19_EP,m19_EP,s19_EP = euler_progresivo(T,dt,alpha19)

#Gráfico
graficarC(t19_EP,s19_EP,"1.9")

#%%

#Lista con los valores de alpha
alphas = [alpha13['alpha'],alpha14['alpha'],alpha15['alpha'],\
          alpha16['alpha'],alpha17['alpha'],alpha18['alpha'],\
          alpha19['alpha']]

#Listas con los valores de los periodos para cada valor en "alphas"
periodos_EP = [periodo(t13_EP,s13_EP),periodo(t14_EP,s14_EP),\
               periodo(t15_EP,s15_EP),periodo(t16_EP,s16_EP),\
               periodo(t17_EP,s17_EP),periodo(t18_EP,s18_EP),\
               periodo(t19_EP,s19_EP)]
    
#%%
    
#Creación del gráfico de "periodos límites en función de alpha"
    
#Creación del lienzo
fig, ax = plt.subplots(figsize=(12,6))

#Gráficos
ax.plot(alphas,periodos_EP)

#Etiquetas
ax.set_xlabel("Valor de $\\alpha$",labelpad=20)
ax.set_ylabel("Periodo límite (millones de años)",labelpad=20)

#Título
ax.set_title("Periodo límite en función de alpha (Método Euler progresivo)",\
             fontweight="bold", loc='center',pad=20)

#Configuraciones
ax.grid(b=True, which='major', axis='both')
ax.margins(0.1)

#Guardado de figura
fig.savefig('Imágenes/C Periodo en función de alpha (EP).pdf')

