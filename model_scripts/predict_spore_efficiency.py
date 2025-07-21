import numpy
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp


# using solve_ivp because its more flexible


r_0 = 10**3
n_v_0 = 10**7
n_s_0 = 0

state_spores_0 = [r_0, n_v_0, n_s_0]
state_0 = [r_0, n_v_0]

r_max = 2
k = 3.2*r_0

s = 10**-7
sigma = 0.282*k
d_max = 1

k_spore = 300.2*r_0
r_max_spore = 1/8

y_v = 10**6 # (maybe)
#(4*10**10)/ (2*10**9) = 20, spore formation is 20X more efficient than cell formation
y_s = 20*y_v

# when the growth cycle ends and CFUs are sampled
t_final = 24

#t_sat = 



def growth_cycle(t, state, y_s):

    #y_s = params[0]
    
    r, n_v = state

    if r < 0:
        r == 0
     
    g_r = r_max * (r/(r+k))

    dr = -1 * (1/y_v)*g_r*n_v 
    dn_v = n_v*g_r 
     
    return [dr, dn_v]


def growth_cycle_spores(state, t, y_s):

    #y_s = params[0]
    
    r, n_v, n_s = state

    if r < 0:
        r == 0
     
    g_r = r_max * (r/(r+k))
    f_r = d_max/ (1 + numpy.exp(s*(r-sigma)))
    g_spore_r = r_max_spore * r/(r+k_spore)
    #g_spore_r = 1

    dr = -1 * (1/y_v)*g_r*n_v  -1 * (1/y_s)*f_r*g_spore_r*n_v
    dn_v = n_v*g_r - n_v*f_r*g_spore_r
    dn_s = 0 + n_v*f_r
     
    return [dr, dn_v, dn_s]




#t = numpy.arange(0.0, t_final, 0.1)
t_span = (0.0, t_final)

y0 = [r_0, n_v_0]
result_ivp = solve_ivp(growth_cycle, t_span, y0, args=(y_s,))

print(result_ivp)



atp_ratio_range = numpy.logspace(-2, 2,  num=50, endpoint=True, base=10.0)

spore_efficiency_all = []

for a in atp_ratio_range:
    
    # 
    y_s = a*y_v

    #params = (y_s)

    #result = odeint(growth_cycle_spores, state_spores_0, t, args=(y_s,))
    #result_ivp = solve_ivp(growth_cycle, t_span, y0, args=(y_s,))

    # find point where resources are depleted 
    #print(sum(result[:,0]<=0))
    #t_idx_final = numpy.argmax(result[:,0]<=0) -1
    #print(t_idx_final)

    #if t_idx_final == 0:
    #    t_idx_final = -1

    #print(r)
    #r_final, n_v_final, n_s_final =  result[t_idx_final,:]
   # r_final, n_v_final, n_s_final =  result[len(t)-1,:]

    #print(min(result[:,0]))

    #spore_efficiency_t = result[:,2] / (result[:,1] + result[:,2])

    #print(r_final)
    #spore_efficiency = n_s_final/(n_s_final+n_v_final)


    #spore_efficiency_all.append(spore_efficiency)



fig, ax = plt.subplots(figsize=(4,4))
ax.plot(atp_ratio_range, spore_efficiency_all)
ax.set_xscale('log', basex=10)
ax.set_xlabel("Ratio of cell vs. spore ATP costs, " + r'$  \frac{\mathrm{ATP}_{\mathrm{cell}} }{ \mathrm{ATP}_{\mathrm{spore}}  }$', fontsize = 9)
ax.set_ylabel("Sporulation efficiency ", fontsize = 9)
ax.axvline(x=20, ls=':', label='Estimate from data', c='k')

ax.legend(loc='upper right', fontsize=7)

#fig.subplots_adjust(hspace=0.4, wspace=0.35)
fig_name = "/Users/williamrshoemaker/GitHub/misc_analyses/test.png"
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()




result_final = odeint(growth_cycle, state_0, t, args=(0.20,))


fig, ax = plt.subplots(figsize=(4,4))

ax.plot(t, result_final[:,0], c='r', label='Resource')
ax.plot(t, result_final[:,1], c='b', label='Vegetative')
#ax.plot(t, result_final[:,2], c='g', label='Spore')

ax.set_yscale('log', basey=10)

print(result_final[:,0])

fig_name = "/Users/williamrshoemaker/GitHub/misc_analyses/trajectory_test.png"
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()


