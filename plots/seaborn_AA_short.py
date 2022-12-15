import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

#labels = ['JC+I+G', 'K2P+I+G', 'F81+I+G', 'F84+I+G','HKY+I+G', 'GTR+I+G']

labels_aa = ['JTT+I+G', 'LG+I+G', 'WAG+I+G', 'Dayhoff+I+G']

sns.set_style("whitegrid")

# Amino acids:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b3,cud,b3,cu)
data = [[37.9,47.43,46.73,47.33],
        [46.23,60.17,60.07,59.93],
        [39.0,38.83,40.03,40.2]]
# Amino acids:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b2,b2,b3,b2)
data_10000 = [[61.17,62.07,59.97,60.63],
        [87.73,87.33,87.0,86.6],
        [50.4,49.07,50.0,49.4]]

x = np.arange(len(labels_aa))  # the label locations
# X = np.arange(len(labels_aa))
Y = np.arange(4)
#figs = plt.figure()
#ax = figs.add_axes([0,0,1,1])

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
width = 0.35  # the width of the bars

# Nucleotides. Sequence length 1000
ax1.yaxis.set_major_formatter('{x:1.0f}%')
ax1.bar(x - width/2, data[0], color = 'darkorange', width=0.15, label='BioNJ')
ax1.bar(x , data[1], color = 'mediumblue', width=0.15, label='ML')
ax1.bar(x + width/2, data[2], color = 'maroon', width=0.15, label='NN')
ax1.set_ylim([0, 100])

# Amino acids. Sequence length 1000
ax2.yaxis.set_major_formatter('{x:1.0f}%')
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
         if counter < 4:
                         ax1.annotate('***',
                                      (bar.get_x() + bar.get_width() / 2,
                                       bar.get_height()), ha='center', va='center',
                                      size=10, xytext=(0, 4),
                                      textcoords='offset points')

         elif counter >= 4 and counter < 8:
                ax1.annotate('***',
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
        if counter < 4:
                 ax2.annotate('***',
                         (bar.get_x() + bar.get_width() / 2,
                          bar.get_height()), ha='center', va='center',
                         size=10, xytext=(0, 4),
                         textcoords='offset points')

        elif counter >= 4 and counter < 8:
                     ax2.annotate('***',
                             (bar.get_x() + bar.get_width() / 2,
                              bar.get_height()), ha='center', va='center',
                             size=10, xytext=(0, 4),
                             textcoords='offset points')

        counter += 1


#fig.suptitle('Accuracy of reconstructions')
#fig.supylabel('Accuracy %')
#fig.legend(shadow=True,fancybox=True)

ax1.set_ylabel('1000 AA')
ax1.set_xticks(x, labels_aa)

ax2.set_ylabel('10000 AA')
ax2.set_xticks(x, labels_aa)


fig.tight_layout()

resolution_value = 200
plt.savefig("/home/nkulikov/Pictures/Final_plots/AA_short_high_resol.png", format="png", dpi=resolution_value)

plt.show()