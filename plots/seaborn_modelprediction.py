import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

labels = ['Nucleotide ML','Nucleotide NN','Aminoacid ML','Aminoacid NN']

labels_aa = ['JTT+I+G LG+I+G WAG+I+G Dayhoff+I+G']

sns.set_style("whitegrid")

# Amino acids:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b3,cud,b3,cu)
data = [87.08,
        92.62,
        100.0,
        99.5]
fig, ax = plt.subplots()
# fig = plt.figure()
#ax = fig.add_axes([0,0,1,1])
ax.yaxis.set_major_formatter('{x:1.0f}%')
plt.bar(labels,data, color = 'maroon', width = 0.4)
ax.set_ylim([0, 100])
# Iterating over the bars one-by-one
counter = 0
for bar in ax.patches:
        # Using Matplotlib's annotate function and
        # passing the coordinates where the annotation shall be done
        # x-coordinate: bar.get_x() + bar.get_width() / 2
        # y-coordinate: bar.get_height()
        # free space to be left to make graph pleasing: (0, 8)
        # ha and va stand for the horizontal and vertical alignment
 #       if counter == 0 or counter == 3 or counter == 6 or counter == 9 or counter == 12 or counter == 15:
             if counter == 0 or counter == 2:
                         ax.annotate('***',
                                      (bar.get_x() + bar.get_width() / 2,
                                       bar.get_height()), ha='center', va='center',
                                      size=12, xytext=(0, 4),
                                      textcoords='offset points')
             counter += 1

#plt.title('Accuracy of reconstructions')
plt.ylabel('Accuracy')
#fig.legend(shadow=True,fancybox=True)



fig.tight_layout()


resolution_value = 200
plt.savefig("/home/nkulikov/Pictures/Final_plots/Model_prediction_final.png", format="png", dpi=resolution_value)

plt.show()