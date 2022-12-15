import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

labels = ['JTT+I+G', 'LG+I+G', 'WAG+I+G','Dayhoff+I+G']

# labels_aa = ['JTT+I+G', 'LG+I+G', 'WAG+I+G', 'Dayhoff+I+G','-','-']
sns.set_style("whitegrid")

# Nucleotides:  data[0] - ML, data[1] -
data = [[0.4580,0.4617,0.4710,0.4824],
        [0.3458,0.3544,0.3465,0.3340]]
# Nucleotides:  data[0] - ML, data[1] -
data_10000 = [ [3.5684,7.2532,22.8530,10.8954],
        [0.3469,0.3498,0.3485,0.3518]]

x = np.arange(len(labels))  # the label locations
# X = np.arange(len(labels_aa))
Y = np.arange(4)
#figs = plt.figure()
#ax = figs.add_axes([0,0,1,1])

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
width = 0.35  # the width of the bars

# Nucleotides. Sequence length 1000
ax1.yaxis.set_major_formatter('{x:1.2f},s')
ax1.bar(x - width/3 , data[0], color = 'mediumblue', width=0.2, label='ML')
ax1.bar(x + width/3, data[1], color = 'maroon', width=0.2, label='NN')
ax1.set_ylim([0, 0.5])
# Amino acids. Sequence length 1000
ax2.yaxis.set_major_formatter('{x:1.2f},s')
ax2.bar(x - width/3, data_10000[0], color = 'mediumblue', width = 0.2)
ax2.bar(x + width/3, data_10000[1], color = 'maroon', width = 0.2)
ax2.set_yscale('log')
ax2.set_ylim([0,10**2])
#ax2.set_ylim([0, 0.5])


#fig.suptitle('Computation time')
#fig.supylabel('seconds')
#fig.legend(shadow=True,fancybox=True)

ax1.set_ylabel('1000 AA')
ax1.set_xticks(x, labels)

ax2.set_ylabel('10000 AA')
ax2.set_xticks(x, labels)


fig.tight_layout()


resolution_value = 200
plt.savefig("/home/nkulikov/Pictures/Final_plots/modelpred_time_aa_high_resol.png", format="png", dpi=resolution_value)

plt.show()
