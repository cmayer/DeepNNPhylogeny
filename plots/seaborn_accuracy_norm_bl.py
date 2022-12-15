import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

labels = ['JC+I+G', 'K2P+I+G', 'F81+I+G', 'F84+I+G','HKY+I+G', 'GTR+I+G']

labels_aa = ['JTT+I+G', 'LG+I+G', 'WAG+I+G', 'Dayhoff+I+G','-','-']

#matplotlib.style.use('fivethirtyeight')

sns.set_style("dark")
# Nucleotides:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b10,b1,b1,b10,cu,cu)
data = [[2174/30,2635/30,2209/30,2619/30,2648/30,2396/30],
        [2986/30,2982/30,2988/30,2984/30,2985/30,2985/30],
        [2984/30,2978/30,2988/30,2983/30,2985/30,2983/30]]
# AminoAcids: data_aa[0] - BioNJ, data_aa[1] - ML, data_aa[2] - NN (b3,b3,b3,b3)
data_aa = [[2777/30,2790/30,2747/30,2775/30,1,1],
        [2998/30,2997/30,2999/30,2997/30,1,1],
        [2865/30,2845/30,2825/30,2821/30,1,1]]

data = np.asarray(data)
data_aa = np.asarray(data_aa)

x = np.arange(len(labels))  # the label locations
X = np.arange(len(labels_aa))
Y = np.arange(4)
#figs = plt.figure()
#ax = figs.add_axes([0,0,1,1])

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
width = 0.35  # the width of the bars

# Nucleotides. Sequence length 1000
ax1.bar(x - width/2, data[0], color = 'lawngreen', width=0.15, label='BioNJ')
ax1.bar(x , data[1], color = 'orangered', width=0.15, label='ML')
ax1.bar(x + width/2, data[2], color = 'turquoise', width=0.15, label='NN')

# Amino acids. Sequence length 1000
ax2.bar(x - width/2, data_aa[0], color = 'lawngreen', width = 0.15)
ax2.bar(x , data_aa[1], color = 'orangered', width = 0.15)
ax2.bar(x + width/2, data_aa[2], color = 'turquoise', width = 0.15)



fig.suptitle('Accuracy of reconstructions')
fig.supylabel('Accuracy %')
fig.legend(shadow=True,fancybox=True)

ax1.set_ylabel('DNA')
ax1.set_xticks(x, labels)

ax2.set_ylabel('Amino acids')
ax2.set_xticks(x, labels_aa)
fig.tight_layout()

#sns.set_style("dark")
#data.plot(kind="bar")
#data_aa.plot(kind="bar")

plt.show()
