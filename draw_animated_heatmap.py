import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors
import matplotlib.animation as animation
genes_covid = [
    (0, 265, ['gene=5UTR'], 0),
    (265, 805, ['gene=nsp1'], 0),
    (805, 2719, ['gene=nsp2'], 0),
    (2719, 8554, ['gene=nsp3'], 0),
    (8554, 10054, ['gene=nsp4'], 0),
    (10054, 10972, ['gene=nsp5'], 0),
    (10972, 11842, ['gene=nsp6'], 0),
    (11842, 12091, ['gene=nsp7'], 0),
    (12091, 12685, ['gene=nsp8'], 0),
    (12685, 13024, ['gene=nsp9'], 0),
    (13024, 13441, ['gene=nsp10'], 0),
    (13441, 13468, ['gene=nsp12'], 0),
    # (13467, 16236, ['gene=nsp12B'], 0),
    # (13475, 13503, ['gene=orf1ab_SL1'], 0),
    # (13487, 13542, ['gene=orf1ab_SL2'], 0),
    (16236, 18039, ['gene=nsp13'], 0),
    (18039, 19620, ['gene=nsp14'], 0),
    (19620, 20658, ['gene=nsp15'], 0),
    (20658, 21552, ['gene=nsp16'], 0),
    # (21552, 21562, ['gene=spacer1'], 0),
    (21562, 25384, ['gene=S'], 0),
    # (25384, 25392, ['gene=spacer2'], 0),

    (25392, 26220, ['gene=ORF3a'], 0),
    # (26220, 26244, ['gene=spacer3'], 0),
    (26244, 26472, ['gene=E'], 0),
    # (26472, 26522, ['gene=spacer4'], 0),

    (26522, 27191, ['gene=M'], 0),
    # (27191, 27201, ['gene=spacer5'], 0),
    (27201, 27387, ['gene=ORF6'], 0),
    # (27387, 27393, ['gene=spacer6'], 0),
    (27393, 27759, ['gene=ORF7'], 0),
    # (27759, 27755, ['gene=spacer7'], 0),
    # (27755, 27887, ['gene=ORF7b'], 0),
    # (27887, 27893, ['gene=spacer7'], 0),
    (27893, 28259, ['gene=ORF8'], 0),
    # (28259, 28273, ['gene=spacer7'], 0),
    (28273, 29533, ['gene=N'], 0),
    # (29533, 29557, ['gene=spacer8'], 0),
    (29557, 29674, ['gene=ORF10'], 0),
    # (29608, 29644, ['gene=ORF10_SL1'], 0),
    # (29628, 29657, ['gene=ORF10_SL2'], 0),
    (29674, 29903, ['gene=3UTR'], 0),
    # (29727, 29768, ['gene=STEMLOOP'], 0),
]

# Generate some example data
np.random.seed(0)
data = pd.DataFrame({
    'mut2mut': np.random.randint(1, 6, size=1000),
    'x': np.random.normal(size=1000),
    'y': np.random.normal(size=1000)
})

# Set up the figure and axis
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot()
hb = ax.hist2d([], [], bins=(20, 20), cmap='viridis')
mut2mut = pd.read_csv('/home/daria/Downloads/summary_table_3nt.csv', header=0)
mut2mut = mut2mut.reindex(sorted(mut2mut.columns), axis=1)
values = mut2mut.filter(like='ed')
# # if more than 5 nulls
# num_of_mut = [row > 18 for row in mut2mut.isnull().sum(axis=1).tolist()]
# indices = [i for i, x in enumerate(num_of_mut) if x]
# mut2mut = mut2mut.drop(indices)

value_cols = mut2mut.filter(like='ed').columns.values.tolist()
mut2mut = mut2mut.fillna(1)
for i in value_cols:
    for j in value_cols[1:]:
        mut2mut[j] = mut2mut[j]/mut2mut[i]


# Define the update function for animation
def update(frame):
    ax.clear()  # Clear the current axis
    current_data = mut2mut  # Select data up to the current frame
    hb = ax.hist2d(current_data['x'], current_data['y'], bins=[np.arange(0, 30000, 300), np.arange(0, 30000, 300)],
           norm = mcolors.PowerNorm(1), cmap='magma', weights=mut2mut[frame].ravel())
    ax.set_title(f'Frame: {frame}')  # Update the title
    ax.set(xlim=(0, 29768))
    ax.set_aspect('equal', adjustable='box')
    for i in genes_covid:
        plt.text(i[0], i[0], i[2][0].split("=")[1].split("'")[0], fontsize=10, color="grey")
    return hb

# Create the animation
num_frames = values.columns.values.tolist()
ani = FuncAnimation(fig, update, frames=num_frames, blit=False, interval=1000)  # Adjust interval as needed
# To save the animation using Pillow as a gif
writer = animation.PillowWriter(fps=15,
                                metadata=dict(artist='None'),
                                bitrate=1800)
ani.save('scatter2021_without_indels.gif', writer=writer)

plt.show()
