'''
Created on Sep 2, 2014

Usage:

Notes:
- Note that paths are relative to the current directory.
- Chromosomes range from 1-38 as the canine genome is used. 
- Exclusion of chrX is arbitrary, it is just as I used in my project

Author: @jaquol
'''



'''Load python modules, general plotting parameters and import files'''

# Reads/deals with table-like files
import pandas as pd
# For generating plots
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# Import file with genomic features
features1 = pd.read_csv('examples/features1.bed', sep = '\t', header = None)
features2 = pd.read_csv('examples/features2.bed', sep = '\t', header = None)
features1.columns = ['chr', 'start', 'end']
features2.columns = ['chr', 'start', 'end']

# Import file with chromosome sizes, sort chromosomes and exclude chrX
chrs = pd.read_csv('examples/chromInfo.bed', sep = '\t', header = None)
chrs.columns = ['chr', 'size']
sorted_chrs = pd.DataFrame(['chr' + str(i) for i in range(1, 39)], columns = ['chr'])
chrs = pd.merge(sorted_chrs, chrs, on = 'chr')
chrs['chr_numeric'] = [int(i.split('chr')[1]) for i in list(chrs['chr'])]



'''Plotting'''

# Plotting parameters
### General plotting parameters
font = {'weight' : 'medium',
        'size'   : 20,
        'family' : 'sans-serif',
        'sans-serif' : ['Arial']}

# Plot chromosomes
plt.close('all')
f, ax = plt.subplots(figsize=(15, 10))
rects = plt.barh(chrs['chr_numeric'], chrs['size'] / 1e6, height = 0.75, color = 'gray', alpha = 0.25, align = 'center')
plt.yticks(chrs['chr_numeric'])
plt.xlabel('Position on chromosome (Mbps)', fontsize = 20)
plt.set_ylim = (0, 125)
yticks_pos = [patch.get_xy()[1] for patch in rects]

# Plot duplications in each individual
a = -1
for chrom in chrs['chr']:  

    # Position counter
    a = a + 1
    
    def plot_duplications(df, offset, my_color):
        
        df = df[df['chr'] == chrom]      
        currentAxis = plt.gca()
        
        for start, end in zip(df['start'], df['end']):
            
            currentAxis.add_patch(Rectangle((start / 1e6, yticks_pos[a] + offset), (end - start) / 1e6, 0.375,
            						facecolor = my_color, alpha = 1, edgecolor = 'none'))
        
    plot_duplications(features1, 0, 'red')
    plot_duplications(features2, 0.375, 'blue')

    plt.text(0.85, 0.95, 'Sample1', transform = ax.transAxes, fontsize = 20, horizontalalignment = 'left', color = 'red')
    plt.text(0.85, 0.90, 'Sample2', transform = ax.transAxes, fontsize = 20, horizontalalignment = 'left', color = 'blue')


plt.savefig('Figures/genomewide.map.duplications.jpg', dpi = 500)