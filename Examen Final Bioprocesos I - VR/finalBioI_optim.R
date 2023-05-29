###################################
#                                 #
# SIMULACION FINAL BIOPROCESOS I  #
# USO DE PKG OPTIM. MAX C_Pf      # 
# MANIUPALDO C_B0                 #
#                                 #
###################################


# CARGA DE LIBRERIAS ----
# = = = = = = = = = = = = 

library(deSolve)



# CONDICIONES DE FERMENTACION ----
# = = = = = = = = = = = = = = = = = =

# ___ Tiempo de fermentacion (simulacion)

t0 = 0 #h 
dt = 2 #h
t1 = 100 #h

times <- seq(t0, t1, by = dt)

n_t = length(times)



# ___ Parametros ----

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



# ___ Condiciones iniciales ----

C_X0 = 1.00695  # g/l
C_S10 = 80  # g/l
C_S20 = 0.25  # g/l
C_P10 = 0  # g/l
C_P20 = 0  # g/l
V0 = 1.5  # l



# DECLARACION DE FUNCIONES ----
# = = = = = = = = = = = = = = =


# ___ mu

mu_crec <- function(mu_m1,C_S1,K_S1,mu_m2,C_S2,K_S2,S_m,n) {
  
  mu = ((mu_m1*C_S1)/(K_S1+C_S1)) * ((mu_m2*C_S2)/(K_S2+C_S2)) * ( 1 - ((C_S1/C_S2)/(S_m))^n )
  
  return(mu)
  
}


# ___ Caudal de alimentacion de glucosa

F_1 <- function(t){
  
  F_1 = (t < 10)*(0) + (t >= 10)*(0.005)
  
  return(F_1)
  
}


# ___ Concentracion de glucosa en la alimentacion

C_S01 <- function(t){
  
  C_S01 = (t < 40)*(300) + (t>=40)*(250)
  
  return(C_S01)
  
}


# ___ Caudal de alimentacion de nitrogeno

F_2 <- function(t){
  
  F_2 = (t < 10)*(0) + (10 < t & t < 40)*(0.005) + (t>=40)*(0)
  
  return(F_2)
  
}


# ___ Caudal total 

F_t <- function(F_1,F_2,t){
  
  F_1 = (t < 10)*(0) + (t >= 10)*(0.005)
  
  F_2 = (t < 10)*(0) + (10 <= t & t < 40)*(0.005) + (t>=40)*(0)
  
  F_t = F_1 + F_2
  
  return(F_t)
  
}



# FUNCION DE RESOLUCION ----
# = = = = = = = = = = = = = =

Bces <- function(t, vals0, parametros) {
  with(as.list(c(vals0, parametros)), {
    
    
    mu = mu_crec(mu_m1,C_S1,K_S1,mu_m2,C_S2,K_S2,S_m,n)
    
    F_1 = F_1(t)
    
    F_2 = F_2(t)
    
    F_t = F_t(F_1,F_2,t)
    
    C_S01 = C_S01(t)
    
    # ___ Balances ----
    # -----------------
    
    dC_X.dt = mu*C_X - ( (F_t/V) * C_X )
    
    dC_S1.dt = -(((mu*C_X)/(Y_XS1)) + m_S1*C_X) + ( ((F_1*C_S01)/V) - ((F_t*C_S1)/V) )
    
    dC_S2.dt = -(((mu*C_X)/(Y_XS2)) + m_S2*C_X) + ( ((F_2*C_S02)/V) - ((F_t*C_S2)/V) )
    
    dC_P1.dt = (K_1*mu*C_X + K_2*C_X) - ((F_t/V)*C_P1)
    
    dC_P2.dt = (K_3*mu*C_X + K_4*C_X) - ((F_t/V)*C_P2)
    
    dV.dt = F_t
    
    
    list(c(dC_X.dt,dC_S1.dt,dC_S2.dt,dC_P1.dt,dC_P2.dt,dV.dt))
    
  })
  
}



# SOLVE ----
# = = = = = =

parametros <- c(mu_m1,K_S1,mu_m2,K_S2,Y_XS1,m_S1,Y_XS2,m_S2,K_1,K_2,K_3,K_4,S_m,n,F_1,F_2,F_t,C_S01,C_S02)
vals0 <- c(C_X = C_X0, C_S1 = C_S10, C_S2 = C_S20, C_P1 = C_P10, C_P2 = C_P20, V = V0)
out <- ode(y=vals0, times = times, func = Bces, parms = parametros)



# POST-PROCESSING ----
# = = = = = = = = = = 

t = out[,1]
res.C_X = out[,2]
res.C_S1 = out[,3]
res.C_S2 = out[,4]*100
#res.C_P1 = out[,5]*10
res.C_P2 = out[,6]*10
res.V = out[,7]



# EJEMPLO OPTIM ----
# = = = = = = = = = =

# ___ Funcion objetivo


f_O <- function(C_X0_p){
  
  vals0 <- c(C_X = C_X0_p, C_S1 = C_S10, C_S2 = C_S20, C_P1 = C_P10, C_P2 = C_P20, V = V0)

  out <- ode(y=vals0, times = times, func = Bces, parms = parametros)

  res.C_P2 = out[,6]*10

  C_P_final = - res.C_P2[length(t)]

  return(C_P_final)
  
}


# ___ Optimizacion

C_X0_optim <- optim(par = C_X0, fn = f_O, method = "L-BFGS-B", 
                                 lower = 0.5*C_X0, upper = 1.5*C_X0)
