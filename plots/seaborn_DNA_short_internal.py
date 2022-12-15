import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

labels = ['JC+I+G', 'K2P+I+G', 'F81+I+G', 'F84+I+G','HKY+I+G', 'GTR+I+G']

# labels_aa = ['JTT+I+G', 'LG+I+G', 'WAG+I+G', 'Dayhoff+I+G','-','-']

sns.set_style("whitegrid")

# Nucleotides:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b10,b1,b1,cud,cu,cu)
data = [[38.43,42.7,39.23,43.3,43.27,43.27],
        [48.0,48.0,47.77,48.33,47.93,48.37],
        [47.57,48.3,48.23,47.93,48.07,47.93]]
# Nucleotides:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b10,b1,b1,cud,cu,cu)
data_10000 = [[40.4,48.97,40.07,49.27,49.17,43.17],
        [73.4,72.87,73.33,73.6,73.17,71.83],
        [71.17,70.23,70.1,71.27,71.0,70.6]]

x = np.arange(len(labels))  # the label locations
# X = np.arange(len(labels_aa))
Y = np.arange(4)
#figs = plt.figure()
#ax = figs.add_axes([0,0,1,1])

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
width = 0.35  # the width of the bars
ax1.yaxis.set_major_formatter('{x:1.0f}%')
# Nucleotides. Sequence length 1000
ax1.bar(x - width/2, data[0], color = 'darkorange', width=0.15, label='BioNJ')
ax1.bar(x , data[1], color = 'mediumblue', width=0.15, label='ML')
ax1.bar(x + width/2, data[2], color = 'maroon', width=0.15, label='NN')
ax1.set_ylim([0, 100])



ax2.yaxis.set_major_formatter('{x:1.0f}%')
# Amino acids. Sequence length 1000

ax2.bar(x - width/2, data_10000[0], color = 'darkorange', width = 0.15)
ax2.bar(x , data_10000[1], color = 'mediumblue', width = 0.15)
ax2.bar(x + width/2, data_10000[2], color = 'maroon', width = 0.15)
ax2.set_ylim([0, 100])


# Iterating over the bars one-by-one
counter = 0
for bar in ax1.patches:
        # Using Matplotlib's annotate function and
        # passing the coordinates where the annotation shall be done
        # x-coordinate: bar.get_x() + bar.get_width() / 2
        # y-coordinate: bar.get_height()
        # free space to be left to make graph pleasing: (0, 8)
        # ha and va stand for the horizontal and vertical alignment
 #       if counter == 0 or counter == 3 or counter == 6 or counter == 9 or counter == 12 or counter == 15:
         if counter < 6:
            ax1.annotate('***',
                       (bar.get_x() + bar.get_width() / 2,
                        bar.get_height()), ha='center', va='center',
                       size=10, xytext=(0, 4),
                       textcoords='offset points')
         elif counter >= 6 and counter < 12:
            ax1.annotate('ns',
                       (bar.get_x() + bar.get_width() / 2,
                        bar.get_height()), ha='center', va='center',
                       size=10, xytext=(0, 4),
                       textcoords='offset points')
         counter += 1
counter = 0

for bar in ax2.patches:
        # Using Matplotlib's annotate function and
        # passing the coordinates where the annotation shall be done
        # x-coordinate: bar.get_x() + bar.get_width() / 2
        # y-coordinate: bar.get_height()
        # free space to be left to make graph pleasing: (0, 8)
        # ha and va stand for the horizontal and vertical alignment
 #       if counter == 0 or counter == 3 or counter == 6 or counter == 9 or counter == 12 or counter == 15:
         if counter < 6:
            ax2.annotate('***',
                       (bar.get_x() + bar.get_width() / 2,
                        bar.get_height()), ha='center', va='center',
                       size=10, xytext=(0, 4),
                       textcoords='offset points')
         elif counter >= 6 and counter < 12:
            if counter == 7 or counter == 9:
                ax2.annotate('*',
                       (bar.get_x() + bar.get_width() / 2,
                        bar.get_height()), ha='center', va='center',
                       size=10, xytext=(0, 4),
                       textcoords='offset points')
            elif counter == 8:
                    ax2.annotate('**',
                                 (bar.get_x() + bar.get_width() / 2,
                                  bar.get_height()), ha='center', va='center',
                                 size=10, xytext=(0, 4),
                                 textcoords='offset points')
            else:
                    ax2.annotate('ns',
                                 (bar.get_x() + bar.get_width() / 2,
                                  bar.get_height()), ha='center', va='center',
                                 size=10, xytext=(0, 4),
                                 textcoords='offset points')
         counter += 1

ax1.set_ylabel('1000 bp')
ax1.set_xticks(x, labels)

ax2.set_ylabel('10000 bp')
ax2.set_xticks(x, labels)


fig.tight_layout()
resolution_value = 200
plt.savefig("myImage.png", format="png", dpi=resolution_value)

plt.show()
#resolution_value = 1200

#plt.savefig("myImage.png", format="png", dpi=resolution_value)
