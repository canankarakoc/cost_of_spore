import numpy
import config

import matplotlib.pyplot as plt


def calculate_deletion_rate():

    indel_summary_file_path = "%sma_summary.csv" % config.data_directory
    indel_summary_file = open(indel_summary_file_path, 'r')

    header_1 = indel_summary_file.readline()
    header_2 = indel_summary_file.readline()
    header_3 = indel_summary_file.readline()

    del_rate_all = []

    for line in indel_summary_file:

        line_split = line.strip().split(',')

        n_del = int(line_split[8])
        n_ins = int(line_split[7])

        n_indels = n_del+n_ins
        if n_indels ==0:
            continue

        rel_del = n_del/(n_del+n_ins)
        indel_rate = float(line_split[-1])
        del_rate_all.append(indel_rate*rel_del)
        #print(line_split)


    indel_summary_file.close()



    print("Mean deletion rate (* 10^-10)/events/gen = %0.2f" % numpy.mean(del_rate_all))



indel_size_file_path = "%sma_size.csv" % config.data_directory
indel_size_file = open(indel_size_file_path, 'r')

del_size_all = []
for line in indel_size_file:

    line_split = line.strip().split(',')

    indel_size = int(line_split[2])

    if indel_size < 0:
        del_size_all.append(abs(indel_size))



del_size_all = numpy.asarray(del_size_all)


# plot CDF

fig, ax = plt.subplots(figsize=(4,4))

size_range = numpy.logspace(0, numpy.log10(max(del_size_all)), num=1000, endpoint=True)
survival_array = numpy.asarray([sum(del_size_all>=s) for s in size_range])/len(del_size_all)
ax.plot(size_range, survival_array, ls='-', lw=1.5, alpha=1, c='k')


ax.set_ylim([0, 1])
ax.set_xscale('log', basex=10)
#ax.set_yscale('log', basey=10)
ax.tick_params(axis='x', labelsize=7)
ax.tick_params(axis='y', labelsize=7)




ax.set_xlabel("Deletion size, " + r'$\Delta$', fontsize = 12 )
ax.set_ylabel("Fraction " + r'$\geq \Delta$', fontsize = 12)



fig.subplots_adjust(hspace=0.60,wspace=0.45)
fig_name = "%sdeletion_size_dist.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()


print(del_size_all)




# plot evo figure
#prob_del_size = numpy.asarray([sum(del_size_all==i)/sum(del_size_all>0) for i in numpy.unique(del_size_all) if (i>0)])





