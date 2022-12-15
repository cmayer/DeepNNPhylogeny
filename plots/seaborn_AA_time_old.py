import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

#labels = ['JC+I+G', 'K2P+I+G', 'F81+I+G','HKY+I+G', 'GTR+I+G']

labels_aa = ['JTT+I+G', 'LG+I+G', 'WAG+I+G', 'Dayhoff+I+G']

sns.set_style("whitegrid")
# WAG: 0.051, 7.9988, 0.383
data = [[0.0493,0.0484,0.051,0.0511],
        [0.7005,0.6901,0.7328,0.8301],
        [0.3849,0.3858,0.383, 0.2701]]
# WAG: 1.6367, 139.6947, 0.3872
data_10000 = [[0.3967,0.3645,1.6367,0.5214],
        [5.4272,13.5457,13.96947,7.624],
        [0.3823,0.3857,0.3872,0.3865]]

x = np.arange(len(labels_aa))  # the label locations
# X = np.arange(len(labels_aa))
Y = np.arange(4)
#figs = plt.figure()
#ax = figs.add_axes([0,0,1,1])

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
width = 0.35  # the width of the bars

# Nucleotides. Sequence length 1000
ax1.yaxis.set_major_formatter('{x:1.2f},s')
ax1.bar(x - width/2, data[0], color = 'darkorange', width=0.15, label='BioNJ')
ax1.bar(x , data[1], color = 'mediumblue', width=0.15, label='ML')
ax1.bar(x + width/2, data[2], color = 'maroon', width=0.15, label='NN')
ax1.set_ylim([0, 1])

# Amino acids. Sequence length 1000
ax2.yaxis.set_major_formatter('{x:1.2f},s')
ax2.bar(x - width/2, data_10000[0], color = 'darkorange', width = 0.15)
ax2.bar(x , data_10000[1], color = 'mediumblue', width = 0.15)
ax2.bar(x + width/2, data_10000[2], color = 'maroon', width = 0.15)
#ax2.set_yscale('log')
#ax2.set_ylim([0,10**2])
#ax2.set_ylim([0, 1])


ax1.set_ylabel('1000 AA')
ax1.set_xticks(x, labels_aa)

ax2.set_ylabel('10000 AA')
ax2.set_xticks(x, labels_aa)


fig.tight_layout()

resolution_value = 200
plt.savefig("/home/nkulikov/Pictures/Final_plots/top_time_AA_high_resol.png", format="png", dpi=resolution_value)

plt.show()
