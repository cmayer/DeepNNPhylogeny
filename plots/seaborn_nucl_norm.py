import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns



labels = ['JC+I+G', 'K2P+I+G', 'F81+I+G', 'F84+I+G','HKY+I+G', 'GTR+I+G']

# labels_aa = ['JTT+I+G', 'LG+I+G', 'WAG+I+G', 'Dayhoff+I+G','-','-']


#matplotlib.style.use('fivethirtyeight')

#sns.set_style("dark")
sns.set_style("whitegrid")
#sns.set_style("white")
#sns.set_style("darkgrid")
#sns.set_style("ticks")

# Nucleotides:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b10,b1,b1,cud,cu,cu)
data = [[72.47,87.83,73.63,87.3,88.27,79.87],
        [99.5,99.4,99.6,99.47,99.5,99.5],
        [99.5,99.27,99.6,99.57,99.5,99.43]]
# Nucleotides:  data[0] - BioNJ, data[1] - ML, data[2] - NN (b10,b1,b1,cud,cu,cu)
data_10000 = [[72.9,87.0,73.33,86.8,86.83,81.77],
        [100.0,100.0,100.0,100.0,100.0,100.0],
        [100.0,100.0,100.0,100.0,99.97,100.0]]

x = np.arange(len(labels))  # the label locations
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

# Amino acids. Sequence length 10000

ax2.yaxis.set_major_formatter('{x:1.0f}%')
ax2.bar(x - width/2, data_10000[0], color = 'darkorange', width = 0.15)
ax2.bar(x , data_10000[1], color = 'mediumblue', width = 0.15)
ax2.bar(x + width/2, data_10000[2], color = 'maroon', width = 0.15)
ax2.set_ylim([0, 100])

#ax1.annotate(text="p-value",xy=(1,120),annotation_clip=False)
#ax1.annotate(text="p-value",xy=(2,100),annotation_clip=False)

#fig.suptitle('Accuracy of reconstructions')
#fig.supylabel('%')

#legend = fig.legend(loc='upper right', frameon= False, ncol=3)
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
plt.savefig("/home/nkulikov/Pictures/Final_plots/DNA_norm_high_resol.png", format="png", dpi=resolution_value)

plt.show()

legend = fig.legend(loc='upper right', frameon= False, ncol=3)
# Save the figure, ensuring that the legend doesn't get cut off
# using `bbox_extra_artists`
figpath = 'graph.png'
fig.savefig(figpath, bbox_extra_artists=[legend], bbox_inches='tight')

# ---------------------------------------------------------------------
# Create a separate figure for the legend
# ---------------------------------------------------------------------
# Bounding box for legend (in order to get proportions right)
# Issuing a draw command seems necessary in order to ensure the layout
# happens correctly; I was getting weird bounding boxes without it.
fig.canvas.draw()
# This gives pixel coordinates for the legend bbox (although perhaps
# if you are using a different renderer you might get something else).
legend_bbox = legend.get_tightbbox(fig.canvas.get_renderer())
# Convert pixel coordinates to inches
legend_bbox = legend_bbox.transformed(fig.dpi_scale_trans.inverted())

# Create teh separate figure, with appropriate width and height
# Create a separate figure for the legend
legend_fig, legend_ax = plt.subplots(figsize=(legend_bbox.width, legend_bbox.height))

# Recreate the legend on the separate figure/axis
legend_squared = legend_ax.legend(
    *ax1.get_legend_handles_labels(),
    bbox_to_anchor=(0, 0, 1, 1),
    bbox_transform=legend_fig.transFigure,
    frameon=False,
    fancybox=None,
    shadow=False,
    ncol=3,
    mode='expand',
)

# Remove everything else from the legend's figure
legend_ax.axis('off')

# Save the legend as a separate figure
legend_figpath = '/home/nkulikov/Pictures/Final_plots/graph-legend_final.png'
print(f"Saving to: {legend_figpath}")
legend_fig.savefig(
    legend_figpath,
    bbox_inches='tight',
    bbox_extra_artists=[legend_squared],
    dpi=resolution_value
)