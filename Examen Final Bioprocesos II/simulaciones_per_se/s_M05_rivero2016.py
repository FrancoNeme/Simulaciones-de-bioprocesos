#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:44:15 2023

########################    
#                      #
# FUNCIONES DEL MODELO #
# M05-rivero2016       #
#                      #
########################

@author: Franco Neme
"""

#%% CARGA DE LIBRERIAS
# = = = = = = = = = = =

import os
import numpy as np
from scipy.integrate import odeint
# import matplotlib.pyplot as plt
import pandas as pd
from graficador_generico import Graficas # Este modulo ya esta adherido al path



#%% FUNCIONES (expresiones intermedias)
# = = = = = = = = = = = = = = = = = = =

def f_mu_x_R(mu_m, N, k_N, k_IN, Gly, k_IG, k_G):
    
    mu_x_R = ( (mu_m * N)/(k_N + N) ) * ( 1 / ( 1 + (N/k_IN) ) ) * ( Gly / (Gly *  (1 + (Gly/k_IG) ) + k_G) ) 
    
    return mu_x_R


def f_mu_p(beta, Gly, k_p, N, k_INP):#, n, PHB, k_d):
    
    mu_p = ( ( beta * Gly) / (Gly + k_p + ( (N)**2 / k_INP) ) ) # * ( 1 - ( (PHB**n) / (PHB**n + k_d)  ) )
    
    return mu_p



#%% CARPETA DE SALIDAS
# = = = = = = = = = =

modelo = 'M05_rivero2016'

Carpeta_salidas = os.path.join(os.getcwd(), '..', 'Salidas', modelo)

try:
    os.mkdir(Carpeta_salidas) # Creo directorio para las salidas del estimador
except:
    pass


#%% SETEO DE LAS CONDICIONES EXPERIMENTALES
# = = = = = = = = = = = = = = = = = = = = =

# Creo clase para importar atributos después

# Tiempo de simulacion
# --------------------

t = np.arange(0,130,1)


# Cond. Inic.
# -----------

x_R0 = 0.25 # g/L
Gly0 = 25 # g/L 
PHB0 = 0 # g/L
N0 = 0.2 # g/L

Balances0 = [x_R0, Gly0, PHB0, N0]


# Cond. Contorno.
# ---------------



#%% PARAMETROS
# = = = = = = =

mu_m = 0.35 # h**−1
k_N = 0.01 # g/L
k_IN = 0.27 # g/L
k_IG = 86.31 # g/L
k_G = 1.05 # g/L
k_x = 9.00 # g/L
n = 1.00 # 
beta = 0.14 # h**−1
k_p = 0.01 # g/L
k_INP = 1E-3 # (g/L)**2
Y_xRGly = 0.35 # g/g
Y_PHBGly = 0.29 # g/g
Y_xRN = 6.0 # g/g
k_d = 9.00 # g/L

parms = mu_m, k_N, k_IN, k_IG, k_G, beta, k_p, k_INP, k_x, Y_xRGly, Y_PHBGly, Y_xRN



#%% RESOLUCION DE SIST. DE ECUACIONES DIF.
# = = = = = = = = = = = = = = = = = = = = =

def f(Balances, mu_m, k_N, k_IN, k_IG, k_G, k_x, beta, k_p, k_INP, Y_xRGly, Y_PHBGly, Y_xRN, t):
    
    # Tupla de salida 
    # ---------------
    
    x_R, Gly, PHB, N  = Balances
    
    # Asignacion de variables acorde a sus funciones 
    # ----------------------------------------------
    
    mu_x_R = f_mu_x_R(mu_m, N, k_N, k_IN, Gly, k_IG, k_G)
    mu_p = f_mu_p(beta, Gly, k_p, N, k_INP)
    
    # Sistema ODE
    # -----------
    
    dx_R_dt = mu_x_R * x_R * (1 - ((x_R + PHB)/k_x)  )
    dGly_dt = - ( (1/Y_xRGly) * mu_x_R + (1/Y_PHBGly) * mu_p ) * x_R * (1 - ((x_R + PHB)/k_x)  )
    dPHB_dt = mu_p * x_R * (1 - ((x_R + PHB)/k_x)  )
    dN_dt = - (1/Y_xRN) * mu_x_R * x_R * (1 - ((x_R + PHB)/k_x)  )
    
    # Salida de la func.
    # ------------------
    
    salida_func =  np.array([dx_R_dt, dGly_dt, dPHB_dt , dN_dt])
    
    return salida_func



#%% ODEINT 
# = = = = =

resultados = odeint(f, Balances0, t, args=parms, mxstep = 5000)


#%% POST-PROCESSING
# = = = = = = = = =

# Variables de estado
# -------------------

df_salidas = pd.DataFrame({'t': t, 'x_R': resultados[:,0], 'Gly': resultados[:,1], 'PHB': resultados[:,2], 'N': resultados[:,3]})

# Funciones internas (tazas de reaccion y x_T)
# --------------------------------------------

x_T = df_salidas['x_R'] + df_salidas['PHB']
mu_x_R = f_mu_x_R(mu_m, df_salidas['N'], k_N, k_IN, df_salidas['Gly'], k_IG, k_G)
mu_p = f_mu_p(beta, df_salidas['Gly'], k_p, df_salidas['N'], k_INP)


# Dataframe de salida general
# ----------------------------

df_salidas = df_salidas.assign(x_T = x_T,
                               mu_x_R = mu_x_R,
                               mu_p = mu_p)


name_df = str(Carpeta_salidas + '/' + 'resultados_simulacion.csv')

df_salidas.to_csv(name_df, index = False, header=True)


#%% Graficas de simulacion
# = = = = = = = = = = = = =

g = Graficas()

for i in df_salidas:
    
    g.grf(t, df_salidas[i], tit_x = 't [h]', tit_y = i + ' [g/L]', fig_name = Carpeta_salidas + '/' + 'sim_' + i)
    

