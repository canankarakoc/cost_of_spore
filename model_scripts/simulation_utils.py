import numpy
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp



colors_dict = {0:'#87CEEB', 1: '#FFA500', 2:'#FF6347'}


# do not change these parameters
n_s_0 = 0
d_max = 1
r_max_spore = 1/8

y_v = 10**6 # (maybe)
#(4*10**10)/ (2*10**9) = 20, spore formation is 20X more efficient than cell formation
#y_s = 20*y_v
y_s = 20*y_v



#r_0 = 10**3
#n_v_0 = 10**7


#state_spores_0 = [r_0, n_v_0, n_s_0]
#state_0 = [r_0, n_v_0]

#r_max = 2
#k = 3.2*r_0

#s = 10**-7
#sigma = 0.282*k


#k_spore = 300.2*r_0
#r_max_spore = 1/8

#y_v = 10**6 # (maybe)
#(4*10**10)/ (2*10**9) = 20, spore formation is 20X more efficient than cell formation
#y_s = 20*y_v

# when the growth cycle ends and CFUs are sampled
t_final = 24
t_span = (0.0, t_final)



def growth_cycle(t, state, y_v, r_max, k):

    #y_s = params[0]
    
    r, n_v = state

    #if r < 0:
    #    r == 0
     
    g_r = r_max * (r/(r+k))

    dr = -1 * (1/y_v)*g_r*n_v 
    dn_v = n_v*g_r 
     
    return [dr, dn_v]



def growth_cycle_spores(t, state, y_v, y_s, r_max, k, d_max, s, sigma, r_max_spore, k_spore):

    #y_s = params[0]
    
    r, n_v, n_s = state

    #if r < 0:
    #    r == 0
     
    g_r = r_max * (r/(r+k))
    f_r = d_max/ (1 + numpy.exp(s*(r-sigma)))
    g_spore_r = r_max_spore * r/(r+k_spore)
    #g_spore_r = 1

    dr = -1 * (1/y_v)*g_r*n_v  -1 * (1/y_s)*f_r*g_spore_r*n_v
    #dn_v = n_v*g_r - n_v*f_r*g_spore_r
    dn_v = n_v*(g_r - f_r*g_spore_r)
    dn_s = 0 + n_v*f_r*g_spore_r
     
    return [dr, dn_v, dn_s]



#result_ivp = solve_ivp(growth_cycle, t_span, state_0, args=(y_s,))
#result_spore_ivp = solve_ivp(growth_cycle_spores, t_span, state_spores_0, args=(y_s,))

