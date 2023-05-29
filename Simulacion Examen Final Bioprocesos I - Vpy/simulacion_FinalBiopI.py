import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd

#%% SETEO DEL SISTEMA

## Parametros ----

mu_m1 = 1.848  # 1/h
K_S1 = 101.78  # g/l
mu_m2 = 0.502  # 1/h
K_S2 = 0.445  # g/l
Y_XS1 = 0.45998  # g/g
m_S1 = 0.0314  # g/(g.h)
Y_XS2 = 44.444  # g/g
m_S2 = 0.00000109  # g/(g.h)
K_1 = 0.00001075  # g/g
K_2 = 0.00000369  # g/(g.h)
K_3 = 0.0000000246  # g/g
K_4 = 0.000520  # g/(g.h)
S_m = 917.8
n = 2.011 
C_S02 = 5.3  # g/l

parms = mu_m1, K_S1, mu_m2, K_S2, Y_XS1, m_S1, Y_XS2, m_S2, K_1, K_2, K_3, K_4, S_m, n, C_S02

## Condiciones iniciales ----

C_X0 = 1.00695  # g/l
C_S10 = 80  # g/l
C_S20 = 0.25  # g/l
C_P10 = 0  # g/l
C_P20 = 0  # g/l
V0 = 1.5  # l

Balances0 = [C_X0, C_S10, C_S20, C_P10, C_P20, V0]

t = list(range(0,102,2));t = np.array(t) 



#%% PLANTEO E INTEGRACION NUMERICA DEL SET DE ECUACIONES DIFERENCIALES

def f(Balances, t, mu_m1, K_S1, mu_m2, K_S2, Y_XS1, m_S1, Y_XS2, m_S2, K_1, K_2, K_3, K_4, S_m, n, C_S02):
    
    C_X, C_S1, C_S2, C_P1, C_P2, V = Balances
    
    # mu

    mu = ((mu_m1*C_S1)/(K_S1+C_S1)) * ((mu_m2*C_S2)/(K_S2+C_S2)) * ( 1 - ((C_S1/C_S2)/(S_m))**n )

    # Caudal de alimentacion de glucosa

    F_1 = np.piecewise(t, [t<10, t>=10], [0, 0.005])

    # Concentracion de glucosa en la alimentacion
    
    C_S01 = np.piecewise(t, [t<40, t>=40], [300, 250])

    # Caudal de alimentacion de nitrogeno

    F_2 = np.piecewise(t, [t<10, (t>=10) & (t < 40), t>=40], [0, 0.005, 0])

    # Caudal total

    F_t = F_1 + F_2

    
    return [mu*C_X - ( (F_t/V) * C_X ), 
           -(((mu*C_X)/(Y_XS1)) + m_S1*C_X) + ( ((F_1*C_S01)/V) - ((F_t*C_S1)/V) ),
           -(((mu*C_X)/(Y_XS2)) + m_S2*C_X) + ( ((F_2*C_S02)/V) - ((F_t*C_S2)/V) ),
           (K_1*mu*C_X + K_2*C_X) - ((F_t/V)*C_P1),
           (K_3*mu*C_X + K_4*C_X) - ((F_t/V)*C_P2),
           F_t]


columnas_salidas = 'Biomasa', 'Glucosa', 'Nitrogeno','Bicaverina','GA3','Volumen' # Me sirve para recordar el orden de salida

# Resuelvo los balances ----

resultados = odeint(f, Balances0, t, args=parms)



#%% POST PROCESSING

# Datos para plot tipo array ----

C_S2P = resultados[:,2]*100
C_P2P = resultados[:,4]*10

resultados_plot1 = np.copy(resultados)
resultados_plot1[:,2] = C_S2P
resultados_plot1[:,4] = C_P2P
resultados_plot1 = np.delete(resultados_plot1, (1,3,5), axis=1)

resultados_plot2 = np.copy(resultados[:,1])


# Conviertiendo datos para plot y salida tipo DataFrame ----


resultados_plot1_df = pd.DataFrame(resultados_plot1);resultados_plot1_df.columns = ['Biomasa','Nitrogeno','GA3']
resultados_plot1_df.insert(0,"Tiempo",t,True)
resultados_plot1_df.insert(4,"Glucosa",resultados_plot2,True)



#%% PLOT OBJETO DEL TIPO DATAFRAME

resultados_plot1_df.plot(x = "Tiempo", y = ["Biomasa","Nitrogeno","GA3"])
ax = resultados_plot1_df.plot(x = "Tiempo", secondary_y=["Glucosa"], mark_right=False)
ax.set_ylabel("Biomasa, Nitrogeno y GA3 [g/L]");
ax.right_ax.set_ylabel("Glucosa [g/L]");
ax.set_xlabel("Tiempo [h]");
ax.get_legend().set_bbox_to_anchor((1.39, 0.7))
ax.get_legend().set_title("Componentes\n")

figure = plt.gcf() # to get the current figure

figure.set_size_inches(8, 6)
figure.savefig("simulacion.png", dpi=100, bbox_inches='tight')



#%% PLOT OBJETO DEL TIPO SALIDA NUMPY

# Create Plot

fig, ax1 = plt.subplots() 
  
ax1.set_xlabel('Tiempo [h]') 
ax1.set_ylabel('Biomasa, GA3 y Nitrogeno [g/L]') 
plot_1 = ax1.plot(t, resultados_plot1) 
ax1.tick_params(axis ='y') 


# Adding Twin Axes

ax2 = ax1.twinx() 
  
ax2.set_ylabel('Glucosa [g/L]') 
plot_2 = ax2.plot(t, resultados_plot2, color = 'magenta') 
ax2.tick_params(axis ='y') 


# Show plot

plt.show()



#%% SALIDA CSV

resultados_df = pd.DataFrame(resultados);resultados_df.columns = columnas_salidas;resultados_df.insert(0,"Tiempo",t,True) 
resultados_df.to_csv('resultados.csv', index=False)