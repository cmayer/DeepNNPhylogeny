import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

labels = ['JC+I+G', 'K2P+I+G', 'F81+I+G','HKY+I+G', 'GTR+I+G']

# labels_aa = ['JTT+I+G', 'LG+I+G', 'WAG+I+G', 'Dayhoff+I+G','-','-']


# Nucleotides:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b10,b1,b1,cud,cu,cu)
data = [[0.0395,0.0385,0.0372,0.0387,0.0482],
        [0.0803,0.081,0.0799,0.0806,0.0986],
        [0.3164,0.3149,0.3136,0.3166,0.3163]]
# Nucleotides:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b10,b1,b1,cud,cu,cu)
data_10000 = [[0.0534,0.0533,0.0495,0.0498,0.0509],
        [0.1219,0.1303,0.1262,0.1304,0.1441],
        [0.3161,0.3117,0.3177,0.318,0.3155]]

data_100000 = [[0.0988,0.1494,0.1001,0.1457,0.1166],
        [0.393,0.5219,0.3948,0.5226,0.5931],
        [0.3234,0.3263,0.324,0.3242,0.3233]]
x = np.arange(len(labels))  # the label locations
# X = np.arange(len(labels_aa))
Y = np.arange(4)
#figs = plt.figure()
#ax = figs.add_axes([0,0,1,1])
sns.set_style("whitegrid")


fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1)
width = 0.35  # the width of the bars

# Nucleotides. Sequence length 1000
ax1.yaxis.set_major_formatter('{x:1.2f},s')
ax1.bar(x - width/2, data[0], color = 'darkorange', width=0.15, label='BioNJ')
ax1.bar(x , data[1], color = 'mediumblue', width=0.15, label='ML')
ax1.bar(x + width/2, data[2], color = 'maroon', width=0.15, label='NN')
ax1.set_ylim([0, 0.6])
# Amino acids. Sequence length 1000
ax2.yaxis.set_major_formatter('{x:1.2f},s')
ax2.bar(x - width/2, data_10000[0], color = 'darkorange', width = 0.15)
ax2.bar(x , data_10000[1], color = 'mediumblue', width = 0.15)
ax2.bar(x + width/2, data_10000[2], color = 'maroon', width = 0.15)
ax2.set_ylim([0, 0.6])

ax3.yaxis.set_major_formatter('{x:1.2f},s')
ax3.bar(x - width/2, data_100000[0], color = 'darkorange', width = 0.15)
ax3.bar(x , data_100000[1], color = 'mediumblue', width = 0.15)
ax3.bar(x + width/2, data_100000[2], color = 'maroon', width = 0.15)
ax3.set_ylim([0, 0.6])

ax1.set_ylabel('1000 bp')
ax1.set_xticks(x, labels)

ax2.set_ylabel('10000 bp')
ax2.set_xticks(x, labels)

ax3.set_ylabel('100000 bp')
ax3.set_xticks(x, labels)

fig.tight_layout()

resolution_value = 200
plt.savefig("/home/nkulikov/Pictures/Final_plots/top_time_DNA_high_resol.png", format="png", dpi=resolution_value)


plt.show()
