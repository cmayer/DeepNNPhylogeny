import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

labels = ['JC+I+G', 'K2P+I+G', 'F81+I+G','HKY+I+G', 'GTR+I+G']

# labels_aa = ['JTT+I+G', 'LG+I+G', 'WAG+I+G', 'Dayhoff+I+G','-','-']
sns.set_style("whitegrid")

# Nucleotides:  data[0] - ML, data[1] - NN (b10,b1,b1,cud,cu,cu)
data = [[0.0390,0.0383,0.0362,0.0373,0.0431],
        [0.3458,0.3472,0.3445,0.3483,0.3427]]
# Nucleotides:  data[0] - ML, data[1] - NN (b10,b1,b1,cud,cu,cu)
data_10000 = [ [0.0481,0.0483,0.0490,0.0493,0.0505],
        [0.3464,0.3450,0.3494,0.3477,0.3494]]

data_100000 = [[0.0990,0.1508,0.0998,0.1508,0.1140],
        [0.3568,0.3533,0.3513,0.3533,0.3538]]

x = np.arange(len(labels))  # the label locations
# X = np.arange(len(labels_aa))
Y = np.arange(4)
#figs = plt.figure()
#ax = figs.add_axes([0,0,1,1])

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1)
width = 0.35  # the width of the bars

# Nucleotides. Sequence length 1000
ax1.yaxis.set_major_formatter('{x:1.2f},s')
ax1.bar(x - width/4 , data[0], color = 'mediumblue', width=0.2, label='ML')
ax1.bar(x + width/3.5, data[1], color = 'maroon', width=0.2, label='NN')
ax1.set_ylim([0, 0.4])
# Amino acids. Sequence length 1000
ax2.yaxis.set_major_formatter('{x:1.2f},s')
ax2.bar(x - width/4, data_10000[0], color = 'mediumblue', width = 0.2)
ax2.bar(x + width/3.5, data_10000[1], color = 'maroon', width = 0.2)
ax2.set_ylim([0, 0.4])

ax3.yaxis.set_major_formatter('{x:1.2f},s')
ax3.bar(x - width/4, data_100000[0], color = 'mediumblue', width = 0.2)
ax3.bar(x + width/3.5, data_100000[1], color = 'maroon', width = 0.2)
ax3.set_ylim([0, 0.4])

#fig.suptitle('Computation time')
#fig.supylabel('seconds')
#fig.legend(shadow=True,fancybox=True)

ax1.set_ylabel('1000 bp')
ax1.set_xticks(x, labels)

ax2.set_ylabel('10000 bp')
ax2.set_xticks(x, labels)

ax3.set_ylabel('100000 bp')
ax3.set_xticks(x, labels)

fig.tight_layout()


resolution_value = 200
plt.savefig("/home/nkulikov/Pictures/Final_plots/modelpred_time_nucl_high_resol.png", format="png", dpi=resolution_value)

plt.show()
