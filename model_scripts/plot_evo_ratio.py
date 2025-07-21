import config
import numpy
import matplotlib.pyplot as plt



# colors_dict = {'0':'#87CEEB', '1': '#FFA500', '2':'#FF6347'}




# mutation estimates obtained from 
# Sung, Way, et al. "Evolution of the insertion-deletion mutation rate across the tree of life." G3: Genes, Genomes, Genetics 6.8 (2016): 2583-2591.
# 10.1534/g3.116.030890

# substitution rate, per-site, per-individual, per-generation
# sem = standard error of the mean
u_sub = 3.41e-10
u_sub_sem = 0.21e-10


# deletion rate, per-event, per-individual, per-generation
u_indel = 1.18e-10
u_idel_sem = 0.09e-10

# in absence of additional information, assume deletion rate is half of the indel rate.
#u_del = 0.5*u_indel
#u_indel = 1.18e-10
u_del = 0.91e-10

Ne_all = [1e6, 1e7, 1e8]

# in unit ATPs
cost_nucleotide = 50
cellular_budet = 2e10
s = cost_nucleotide/cellular_budet

def calculate_ratio_rates(Ne, deletion_size, u_del_):

    scaled_s = 2*Ne*s*deletion_size

    return (scaled_s/(1-numpy.exp(-scaled_s))) * (u_del_/u_sub)



deletion_size_all = numpy.logspace(0, 3, num=100, endpoint=True, base=10.0)
deletion_size_all = numpy.unique(deletion_size_all.astype(int))


Ne_lower = 6.119e7
Ne_upper = 3.224e8

#ratio_rates_6 = calculate_ratio_rates(1e6, deletion_size_all, u_del)
#ratio_rates_7 = calculate_ratio_rates(1e7, deletion_size_all, u_del)
#ratio_rates_8 = calculate_ratio_rates(1e8, deletion_size_all, u_del)

ratio_rates_lower = calculate_ratio_rates(Ne_lower, deletion_size_all, u_del)
ratio_rates_higher = calculate_ratio_rates(Ne_upper, deletion_size_all, u_del)




fig, ax = plt.subplots(figsize=(4,4))
#ax.plot(deletion_size_all, ratio_rates_6, ls='-', lw=2, alpha=1, label=r'$N_{e} = 10^{6}$', c='#FF6347', zorder=2)
#ax.plot(deletion_size_all, ratio_rates_7, ls='-', lw=2, alpha=1, label=r'$N_{e} = 10^{7}$', c='#FFA500', zorder=2)
#ax.plot(deletion_size_all, ratio_rates_8, ls='-', lw=2, alpha=1, label=r'$N_{e} = 10^{8}$', c='#87CEEB', zorder=2)


ax.plot(deletion_size_all, ratio_rates_lower, ls='-', lw=2, alpha=1, label=r'$N_{e} = 6.1 \cdot 10^{7}$' + ', Sung et al. (2016)', c='#FF6347', zorder=2)
ax.plot(deletion_size_all, ratio_rates_higher, ls='-', lw=2, alpha=1, label=r'$N_{e} = 3.2 \cdot 10^{8}$' + ', Bobay and Ochman (2018)', c='#87CEEB', zorder=2)



ax.axhline(y=u_del/u_sub, lw=2, ls=':', c='k', label='Neutrality')


ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)
ax.set_xlabel("Deletion size (bp), " + r'$\Delta$', fontsize = 12)
ax.set_ylabel("Deletion vs. substitution fixation rates under\nrelaxed selection for spore formation, " + r'$ \frac{d_{\mathrm{del}}}{d_{\mathrm{sub}}}$', fontsize = 11)



ax.legend(loc='upper left', fontsize=7)

#fig.subplots_adjust(hspace=0.4, wspace=0.35)
fig_name = "%sevo_ratio.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



#print(u_del/u_sub)

# save data for Canan as flat file
data_file_path = '%sevo_ratio.tsv' % config.data_directory
data_file = open(data_file_path, 'w')

deletion_size_str =  'deletion_size\t' + ",".join(str(x) for x in deletion_size_all)
ratio_rates_lower_ne_str = 'ratio_rates_lower_Ne\t' + ",".join(str(x) for x in ratio_rates_lower)
ratio_rates_higher_ne_str = 'ratio_rates_higher_Ne\t' + ",".join(str(x) for x in ratio_rates_higher)

data_file.write(deletion_size_str)
data_file.write('\n')
data_file.write(ratio_rates_lower_ne_str)
data_file.write('\n')
data_file.write(ratio_rates_higher_ne_str)
data_file.close()





#del_size_dist = [1,0,1,0,3,0,3,0,1,1,0,2,2,3,4,2,3,3,3,3,2,3,3,3,2,2,4,2,1,3,1,6,3,2,1,1,2,1,0,2,2,1,2,2,1,2,2,1,2,2]

del_size_dist = [1,6,66,1,1,1,10,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,1,1,1,1,1,1,12,9,1,1,1,1,5,1,1,1,1,1,309,1,1401,1401,1401,1401,1401,1401,1401,1401,1401,4136,1,1,1,1,1,2,1,1,1,1,1,1,2,1,1,1,2,1,2,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,21,1,44,39,1,1]
del_size_dist = numpy.asarray(del_size_dist)

prob_del_size = [sum(del_size_dist==i)/sum(del_size_dist>0) for i in numpy.unique(del_size_dist) if (i>0)]
prob_del_size = numpy.asarray(prob_del_size)


deletion_size_conditional_all = numpy.unique(del_size_dist)



#ratio_rates_conditional_6 = calculate_ratio_rates(1e6, deletion_size_conditional_all, u_del*prob_del_size)
#ratio_rates_conditional_7 = calculate_ratio_rates(1e7, deletion_size_conditional_all, u_del*prob_del_size)
#ratio_rates_conditional_8 = calculate_ratio_rates(1e8, deletion_size_conditional_all, u_del*prob_del_size)

ratio_rates_conditional_lower = calculate_ratio_rates(Ne_lower, deletion_size_conditional_all, u_del*prob_del_size)
ratio_rates_conditional_higher = calculate_ratio_rates(Ne_upper, deletion_size_conditional_all, u_del*prob_del_size)


fig, ax = plt.subplots(figsize=(4,4))

#ax.plot(deletion_size_conditional_all, ratio_rates_conditional_6, ls='-', lw=2, alpha=1, label=r'$N_{e} = 10^{6}$', c='#FF6347', zorder=2)
#ax.plot(deletion_size_conditional_all, ratio_rates_conditional_7, ls='-', lw=2, alpha=1, label=r'$N_{e} = 10^{7}$', c='#FFA500', zorder=2)
#ax.plot(deletion_size_conditional_all, ratio_rates_conditional_8, ls='-', lw=2, alpha=1, label=r'$N_{e} = 10^{8}$', c='#87CEEB', zorder=2)


#ax.plot(deletion_size_conditional_all, ratio_rates_conditional_8, ls='-', lw=2, alpha=1, label=r'$N_{e} = 10^{8}$', c='#87CEEB', zorder=2)
ax.plot(deletion_size_conditional_all, ratio_rates_conditional_lower, ls='-', lw=2, alpha=1, label=r'$N_{e} = 6.1 \cdot 10^{7}$' + ', Sung et al. (2016)', c='#FF6347', zorder=2)
ax.plot(deletion_size_conditional_all, ratio_rates_conditional_higher, ls='-', lw=2, alpha=1, label=r'$N_{e} = 3.2 \cdot 10^{8}$' + ', Bobay and Ochman (2018)', c='#87CEEB', zorder=2)
ax.plot(deletion_size_conditional_all,( u_del*prob_del_size)/u_sub, ls=':', lw=2, alpha=1, label="Neutrality", c='k', zorder=2)


null_ratio = (u_del*prob_del_size)/u_sub
ax.scatter(deletion_size_conditional_all, ratio_rates_conditional_lower, s=10, alpha=1, c='#FF6347', zorder=3)
ax.scatter(deletion_size_conditional_all, ratio_rates_conditional_higher, s=10, alpha=1, c='#87CEEB', zorder=3)
ax.scatter(deletion_size_conditional_all, null_ratio, s=10, alpha=1, c='k', zorder=3)

#ax.axhline(y=u_del/u_sub, lw=2, ls=':', c='k', label='Neutrality')

#ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)
ax.set_xlabel("Deletion size (bp), " + r'$\Delta$', fontsize = 12)
ax.set_ylabel("Deletion vs. substitution fixation rate nunder\nrelaxed selection for spore formation, " + r'$ \frac{d_{\mathrm{del}} (\Delta) }{d_{\mathrm{sub}}}$', fontsize = 11)
ax.set_xscale('log', basex=10)
ax.legend(loc='upper left', fontsize=7)

#fig.subplots_adjust(hspace=0.4, wspace=0.35)
#fig_name = "/Users/williamrshoemaker/GitHub/sporecosts_model/analysis/evo_ratio_conditional.png"
fig_name = "%sevo_ratio_conditional.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



# save data for Canan as flat file
data_file_path = '%sevo_ratio_conditional.tsv' % config.data_directory
data_file = open(data_file_path, 'w')

deletion_size_str =  'deletion_size\t' + ",".join(str(x) for x in deletion_size_conditional_all)
ratio_rates_lower_ne_str = 'ratio_rates_lower_Ne\t' + ",".join(str(x) for x in ratio_rates_conditional_lower)
ratio_rates_higher_ne_str = 'ratio_rates_higher_Ne\t' + ",".join(str(x) for x in ratio_rates_conditional_higher)
ratio_rates_higher_ne_str = 'ratio_rates_higher_Ne\t' + ",".join(str(x) for x in ratio_rates_conditional_higher)
null_ratio_rates = 'null_ratio_rates\t' + ",".join(str(x) for x in null_ratio.tolist())

data_file.write(deletion_size_str)
data_file.write('\n')
data_file.write(ratio_rates_lower_ne_str)
data_file.write('\n')
data_file.write(ratio_rates_higher_ne_str)
data_file.write('\n')
data_file.write(null_ratio_rates)
data_file.close()




