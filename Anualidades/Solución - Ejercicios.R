
# ========== #
# EJERCICIOS # 
# ========== #

# LECTURA DE ARCHIVO DE TABLA

options(scipen=999)
TablaSBS = read.table("TMSPP2017.txt",T)
attach(TablaSBS)
library(lifecontingencies)

# EJERCICIO 1

Tabla1 = probs2lifetable(probs = SPPS2017H*(1-AaxH)^(2), 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla1")
TEM = 0.6/100
TEA = (1+TEM)^12-1
a34 = axn(Tabla1, x = 34, i = TEA) - 1
12000*a34 # 153030.7

# EJERCICIO 2

Tabla2 = probs2lifetable(probs = SPPI2017H*(1-AaxH)^(2), 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 2")
a19.50 = axn(Tabla2, x = 19, n = 50, i = 0.02)
18000*a19.50

# EJERCICIO 3

Tabla3 = probs2lifetable(probs = SPPS2017M*(1-AaxM)^(2), 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 3")
TEM    = (5/12)/100
TEA    = (1+TEM)^12-1
a      = axn(Tabla3, x = 31, n = 32, i = TEA)
npx    = pxt(Tabla3, x = 31, t = 32)
a3132  = a - 1 + 1/(1+TEA)^32*npx
a3132*15000 # 231484.5

# EJERCICIO 4

Tabla4 = probs2lifetable(probs = SPPS2017H*(1-AaxH)^(2), 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 4")
a = axn(Tabla4, x = 35, m = 30, n = 20, i = 0.01)
a*15000 # 160754.8

a = Exn(Tabla4, x = 35, n = 30, i = 0.01)*axn(Tabla4, x = 65, n = 20, i = 0.01)
a*15000 # 160754.8

a = axn(Tabla4, x = 35, n = 50, i = 0.01) - axn(Tabla4, x = 35, n = 30, i = 0.01)
a*15000

# EJERCICIO 5

Tabla5 = Tabla4
E  = Exn(Tabla5, x = 35, n = 30, i = 0.01)
a  = axn(Tabla5, x = 65, n = 20, i = 0.01)
p  = pxt(Tabla5, x = 65, t = 20)
aa = E*(a-1+1/1.01^20*p)
aa*15000  
  
a1 = axn(Tabla5, x = 35, n = 50, i = 0.01)
p1 = pxt(Tabla5, x = 35, t = 50)
a2 = axn(Tabla5, x = 35, n = 30, i = 0.01)
p2 = pxt(Tabla5, x = 35, t = 30)
aa = a1 + 1/1.01^50*p1 - a2 - 1/1.01^30*p2
aa*15000  

# EJERCICIO 6

Tabla6 = probs2lifetable(probs = SPPI2017M*(1-AaxM)^(2), 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 6")
a = axn(Tabla6, x = 43, m = 5, i = sqrt(1.07)-1) # considerando esquema anticipado
a
E   = Exn(Tabla6, x = 43, n = 4, i = sqrt(1.07)-1) #considerando esquema vencido
a47 = axn(Tabla6, x = 47, i = sqrt(1.07)-1)
E*(a47-1)
a1 = axn(Tabla6, x = 43, i = sqrt(1.07)-1)
a2 = axn(Tabla6, x = 43, n = 4, i = sqrt(1.07)-1)
 p = pxt(Tabla6, x = 43, t = 4)
a1 - a2 - 1/(1+sqrt(1.07)-1)^4*p

# EJERCICIO 7

TEM = (0.5/3)/100
TEA = (1+TEM)^12-1
d   = log(1+TEA)
S = function(x){(1-x^2/100^2)^0.4}  # fn de supervivencia
f = function(t){exp(-d*t)*S(46.5+t)/S(46.5)} # fn a integrar
integrate(f,9.5,35)$value
integrate(f,9.5,35)$value*800*52

# EJERCICIO 8

Tabla8 = Tabla6
library(FinCal)
a1 = pv.annuity(r = 0.01, n = 15, pmt = -1, type = 1)
a2 = axn(Tabla8, x = 35, m = 15, i = 0.01)
a  = a1+a2 # 34.01956
a*15000 # 510293.4

# PC2
# CAP 4 - SEGUROS DE VIDA (DOTAL PURO, VITALICIO, TEMPORAL, DOTAL, DIFERIDO)
# CAP 5 - ANUALIDADES:
# VITALICIAS: VENCIDA, ANTICIPADA, CONTINUA
# TEMPORALES: VENCIDA, ANTICIPADA, CONTINUA
# SOLO INMEDIATAS, NO DIFERIDAS

