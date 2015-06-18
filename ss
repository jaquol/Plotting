
# coding: utf-8

# # Rearrangement breakpoints
# ------
# 
# The T47D genome has suffered many chromosomal rearrangements compared to the normal human genome. I will use the Lumpy tool to identiy them.
# 
# - Download and Install Lumpy
# - Pre-processing
# - Run lumpyexpress
# - Post-processing
# - Analysis

# ## Download and Install Lumpy
# -----
# 
# [Lumpy](https://github.com/arq5x/lumpy-sv#installation) is a general probabilistic framework for structural variant discovery. I download and installed as follows:
# 
# *In the CRG iMac:*
# 
# > ll /Users/jquilez/Downloads  
# > git clone --recursive https://github.com/arq5x/lumpy-sv.git   
# > mv lumpy-sv /Volumes/users-GR-mb-jquilez  
# 
# 
# *In the CRG Cluster:*
# > mv -v lumpy-sv /software/mb/el6.3/  
# > cd /software/mb/el6.3/lumpy-sv/
# 
# > make  
# > cp -R lumpy-sv /Volumes/users-GR-mb-jquilez/  
# 
# 
# > mv lumpy-sv /software/mb/el6.3/  
# > ln -s /software/mb/el6.3/lumpy-sv/bin/lumpy /software/mb/bin/lumpy  
# > ln -s /software/mb/el6.3/lumpy-sv/bin/lumpyexpress /software/mb/bin/lumpyexpress  
# 
# Modify the '/software/mb/el6.3/lumpy-sv/bin/lumpyexpress.config' file so that it contains all the paths needed:
# 
# > ...  
# > SAMTOOLS=/usr/bin/samtools  
# > ...  

# ## Pre-processing
# -----
# 
# In the Lumpy terminology, this step refers to the alignment of reads and the extraction of dicordant paired-end alignmnents and the split-read alignments. Alignment is not needed here as I previously aligned reads with BWA-MEM (see 'mapping-hg19-bwa.ipynb'), the aligner recommended in the Lumpy workflow. I will use the BAM files after merging sequencing lanes and removing duplicates.
# 
# 
# ### Extract discordant paired-end and split-read alignments
# ------
# 
# Run the script:
# 
# > ./Scripts/lumpy.pre.processing.sh
# 
# Check if the output files are sorted by genomic coordinates:
# 
# > cat Data/Lumpy/job_out/preprocessing.803369.out
# 
# 
# ### Define regions with excessive coverage
# ---
# 
# In the [lumpy paper](http://www.genomebiology.com/2014/15/6/R84) regions with excessive coverage are excluded, with these defined as:
# 
# ```
# depth coverage > (2 * mode) + 3 standard deviations
# ```
# 
# Following this criterion, >52x is defined as excessive coverage in our region (see 'mapping-hg19-bwa.ipynb'). I will run lumpy with/without removing regions showing excessive coverage

# ## Run lumpyexpress
# -----
# 
# Run the following scripts (the second excludes regions with excessive coverage, as defined above):
# 
# > ./Scripts/lumpy.run.lumpyexpress.sh
# > ./Scripts/lumpy.run.lumpyexpress.limit.coverage.sh

# ## Post-processing
# ---
# 
# I download [svtyper](https://github.com/cc2qe/svtyper/edit/master/svtyper) and saved as:
# 
# > Scripts/svtyper
# 
# Run the scritp:
# 
# > ./Scripts/lumpy.svtyper.sh
# 

# ## Analysis
# ----

# ### Python libraries and default parameters
# ----

# In[1]:

get_ipython().magic(u'pylab inline')

import os
import pandas as pd
import glob
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter as scf
import random
from matplotlib import gridspec
from matplotlib.patches import Rectangle
from scipy.stats.stats import pearsonr
from scipy import stats
import pybedtools
from pybedtools import BedTool
import matplotlib_venn
from matplotlib_venn import venn2

# FONT
plt.rcParams['font.size'] = 20 
plt.rcParams['font.weight'] = 'medium' 
plt.rcParams['font.family'] = 'sans-serif' 
plt.rcParams['font.sans-serif'] = 'Arial' 

# LINES
plt.rcParams['lines.linewidth'] = 2.0

# LEGEND
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['legend.frameon'] = False

# SAVING FIGURES
#plt.rcParams['savefig.dpi'] = 300 
plt.rcParams['savefig.bbox'] = 'tight'


# In[2]:

# Set working directory
project = 'T47D_gDNA'
DIR = '/Volumes/users-GR-mb-jquilez/%s' % project


# ### Import data
# ---

# In[3]:

def import_bedpe(limit_coverage):
    
    if limit_coverage == 'no_coverage_limit':
        file_name = '3128_T47D.breakpoints.bedpe'
    elif limit_coverage == '52':
        file_name = '3128_T47D.breakpoints.limit.coverage.%s.bedpe' % limit_coverage
    
        
    infile = '%s/Data/Lumpy/run_lumpyexpress/BEDPE/%s' % (DIR, file_name)
    df = pd.read_table(infile, comment = "#", header = None)
    df.columns = ['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B', 'ID', 'QUAL', 'STRAND_A',
                  'STRAND_B', 'TYPE', 'FILTER', 'INFO', 'FORMAT', '3128_T47D']
    return df


# In[4]:

svs = {}
svs['no_coverage_limit'] = import_bedpe('no_coverage_limit')
svs['52'] = import_bedpe('52')


# ### Overlap calls including/excluding regions with excessive coverage
# ---
# 
# ```
# a=Data/Lumpy/run_lumpyexpress/BEDPE/3128_T47D.breakpoints.bedpe
# b=Data/Lumpy/run_lumpyexpress/BEDPE/3128_T47D.breakpoints.limit.coverage.52.bedpe
# 
# # Number of variants including all regions = 35778
# grep -v CHROM $a | wc -l
# 
# # Number of variants in regions with coverage <52x = 10959
# grep -v CHROM $b | wc -l
# 
# # Number of variants in which both ends overlap in both call sets = 10949
# bedtools pairtopair -a $a -b $b | wc -l
# 
# ```

# In[5]:

def overlap_bedpe():
    
    # Values
    a = 35778
    b = 10959
    a_and_b = 10949
     
    # Plot
    plt.close('all')
    fig, ax = plt.subplots(figsize = (8, 8))
    venn2(subsets = (a - a_and_b, b - a_and_b, a_and_b), set_labels = ('All', '<52x'))
    outfile = '%s/Figures/venn.diagram.overlap.exclude.include.regions.high.coverage.png' % DIR
    savefig(outfile, bbox_inches = 'tight', dpi = 400)   
    
    # Print summary statistics
    print "Not excluding regions with coverage >52x identified %f times more sites" % ((a*1.)/b)
    print "Virtually all regions in the more-restrictive call set are included in the more flexible one"


# In[6]:

overlap_bedpe()


# ### Support from discordant paired-ends and split reads
# -----

# In[27]:

def support_pe_sr(limit_coverage, ax1, ax2, ax3):
    
    # Retrieve data
    df = svs[limit_coverage].copy()
    
    # Extract the number of reads supporting from each class (i.e. discordant paired-ends or split reads)
    # supporting each variant
    df['reads_su'] = [int(i.split(':')[1]) for i in df['3128_T47D']]
    df['reads_pe'] = [int(i.split(':')[2]) for i in df['3128_T47D']]
    df['reads_sr'] = [int(i.split(':')[3]) for i in df['3128_T47D']]

    limit_coverage_to_text = {}
    limit_coverage_to_text['no_coverage_limit'] = 'All'
    limit_coverage_to_text['52'] = 'Coverage <52x'
    
    # Histogram
    def plot_hist(support, color):  
        ax1.hist(df[support], bins = arange(0, 100), normed = True, color = color, histtype = 'step',
                 label = support, linewidth = 2, alpha = 0.85, cumulative = True)
    
    plot_hist('reads_pe', 'blue')
    plot_hist('reads_sr', 'orange')
    plot_hist('reads_su', 'green')
    ax1.set_ylim(0, 1)
    ax1.set_xlim(0, 20)
    ax1.set_xlabel('Number of reads')
    ax1.set_ylabel('Fraction of calls')
    ax1.legend(fontsize = 14, loc = 'lower right')
                   
    # Scatter plot
    ax2.scatter(log10(df['reads_su']), df['reads_pe'] / df['reads_su'], color = 'gray', alpha = 0.10)
    ax2.set_ylim(-0.1, 1.1)
    ax2.set_xlim(0, 5)
    ax2.set_xlabel('log10(reads_pe + reads_sr)')
    ax2.set_ylabel('Fraction of reads_pe')
    ax2.text(0.5, 1.2, limit_coverage_to_text[limit_coverage], transform = ax2.transAxes, fontsize = 22,
            horizontalalignment = 'center')        
    
    # Barplot
    values = []
    ind = arange(3)
    width = 0.25
    for r in arange(3):
        
        cond1 = df['reads_pe'] > r
        cond2 = df['reads_sr'] > r   
        values.append((1. * len(df[cond1])) / len(df))
        values.append((1. * len(df[cond2])) / len(df))
        values.append((1. * len(df[cond1 & cond2])) / len(df))
    
    rects1 = ax3.bar(ind, values[0:3], width, color = 'blue')
    rects2 = ax3.bar(ind + width, values[3:6], width, color = 'orange')
    rects3 = ax3.bar(ind + 2*width, values[6:], width, color = 'green')
    ax3.set_ylabel('Fraction of calls')
    ax3.set_xticks(ind + width)
    ax3.set_ylim(0, 1)
    ax3.set_xticklabels( ('>0', '>1', '>2') )
    ax3.set_xlabel('Number of reads\nsupporting the call')
    ax3.legend( (rects1[0], rects2[0], rects3[0]), ('pe', 'sr', 'pe + sr'), fontsize = 14, bbox_to_anchor=(1.5, 1) )


# In[28]:

plt.close('all')
f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(15, 8))
f.subplots_adjust(wspace = 0.50, hspace = 1)
support_pe_sr('no_coverage_limit', ax1, ax2, ax3)
support_pe_sr('52', ax4, ax5, ax6)
outfile = '%s/Figures/structural.variants.calls.suppor.pe.sr.png' % DIR
#savefig(outfile, bbox_inches = 'tight', dpi = 400)   


# ### Types of variants
# ----

# In[405]:

def types_of_variants(limit_coverage, offset, alpha):
    
    # Retrieve data
    df = svs[limit_coverage].copy()
    
    # Barplot
    ind = arange(4)
    width = 0.35
    variant_totals = []
    for variant_type in set(df['TYPE']):
        
        variant_totals.append(len(df[df['TYPE'] == variant_type]))
    
    rel_freqs = [(1.*i)/sum(variant_totals) for i in variant_totals]
    rects = ax1.bar(ind + offset, rel_freqs, width, color = 'gray', alpha = alpha)
    
    # Axes and labels
    ax1.set_ylabel('Fraction of calls')
    ax1.set_xticks(ind + width)
    ax1.set_ylim(0, 1)
    ax1.set_xlim(-0.1, 4.1)
    ax1.set_xticklabels( list(set(df['TYPE'])) )
    ax1.set_xlabel('Type of structural variant')
    
    # More labels
    labels = {}
    labels['DUP'] = 'Duplication'
    labels['INV'] = 'Inversion'
    labels['DEL'] = 'Deletion'
    labels['BND'] = 'Complex Rearrangement'
    
    for a, b in zip(arange(4), [1., 0.9, 0.8, 0.7]):
        
        ax1.text(ind[a] + offset + 0.05, 0.05 + (1. * variant_totals[a]) / sum(variant_totals), variant_totals[a],
                 fontsize = 16, horizontalalignment = 'left')
        
        if limit_coverage == '52':
            ax1.text(4.5, b, '%s = %s' % (list(set(df['TYPE']))[a], labels[list(set(df['TYPE']))[a]]))
        
    return rects


# In[406]:

plt.close('all')
f, (ax1) = plt.subplots(1, 1, figsize=(10, 4))
rects1 = types_of_variants('no_coverage_limit', 0, 0.85)
rects2 = types_of_variants('52', 0.35, 0.25)
ax1.legend((rects1[0], rects2[0]), ('All', '<52x'), fontsize = 18, loc = 'upper left')
outfile = '%s/Figures/variant.types.png' % DIR
savefig(outfile, bbox_inches = 'tight', dpi = 400)  


# ### Confidence interval of breakpoints
# ----

# In[101]:

def confidence_interval_breakpoints(limit_coverage, ax1, ax2, ax3):
    
    # Retrieve data
    df = svs[limit_coverage].copy()
    df['reads_su'] = [int(i.split(':')[1]) for i in df['3128_T47D']]
     
    labels = {}
    labels['DUP'] = 'blue'
    labels['INV'] = 'orange'
    labels['DEL'] = 'green'
    labels['BND'] = 'brown'

    limit_coverage_to_text = {}
    limit_coverage_to_text['no_coverage_limit'] = 'All'
    limit_coverage_to_text['52'] = 'Coverage <52x'
    
    # Calculate confidence interval for both breakends of the structural variants
    df['end_ci'] = [int(i.split(':')[2]) for i in df['3128_T47D']]

    df['a_ci'] = df['END_A'] - df['START_A']
    df['b_ci'] = df['END_B'] - df['START_B']
    
    # Confidence interval distribution
    for variant_type, i in zip(set(df['TYPE']), arange(1, 9, 2)):
        
        k = df[df['TYPE'] == variant_type]
        myhist = ax1.hist(list(k['a_ci']), histtype = 'step', bins = arange(1000),
                          label = variant_type + '-start', normed = True, cumulative = True,
                          alpha = 0.75, linewidth = 2, color = labels[variant_type])
        myhist = ax1.hist(list(k['b_ci']), histtype = 'step', bins = arange(1000),
                          label = variant_type + '-end', normed = True, cumulative = True,
                          alpha = 0.75, linewidth = 2, color = labels[variant_type],
                             linestyle = 'dashed')
        
    ax1.set_xlim(0, 600)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel('Breakpoint CI (bp)')
    ax1.set_ylabel('Fraction of calls')
    
    # CI vs reads support
    ax2.scatter(log10(df['reads_su']), df['a_ci'], color = 'gray', alpha = 0.1)
    ax2.scatter(log10(df['reads_su']), df['b_ci'], color = 'gray', alpha = 0.1)
    ax2.set_xlim(0, 5)
    ax2.set_ylim(0, 700)
    ax2.set_xlabel('log10(reads pe + reads sr)')
    ax2.set_ylabel('Breakpoint CI (bp)')
    ax2.text(0.5, 1.2, limit_coverage_to_text[limit_coverage], transform = ax2.transAxes, fontsize = 22,
            horizontalalignment = 'center')   
    
    # CI vs variant size (for DUP, INV, DEL variants only --BND do not have size)
    df_sub = df[df['TYPE'] != 'BND'].copy()
    df_sub['SVLEN'] = [abs(int(i.split(';')[1].split('=')[1])) for i in df_sub['INFO']]
    ax3.scatter(log10(df_sub['SVLEN']), df_sub['a_ci'], color = 'gray', alpha = 0.1)
    ax3.scatter(log10(df_sub['SVLEN']), df_sub['b_ci'], color = 'gray', alpha = 0.1)
    ax3.set_xlim(0, 10)
    ax3.set_ylim(0, 700)
    ax3.set_xlabel('log10(variant size)')
    ax3.set_ylabel('Breakpoint CI (bp)')

    if limit_coverage == '52':
        ax1.legend(loc = 'lower right', fontsize = 10, ncol = 2)


# In[102]:

plt.close('all')
f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(16, 10))
f.subplots_adjust(wspace = 0.50, hspace = 1)
confidence_interval_breakpoints('no_coverage_limit', ax1, ax2, ax3)
confidence_interval_breakpoints('52', ax4, ax5, ax6)


# ## Summary
# ----
# 
# - Not excluding regions with coverage >52x identified 3.264714 times more sites; virtually all regions in the more-restrictive call set (coverage<52x) are included in the more flexible one
# - ~70% of variant calls are supported by 2 discordant paired-ends (pe) or 2 split reads (sr), while the rest of the variants are supported by a broads range of reads
# - More than 80% of variants are supported by either pe or sr, with <25% supported by both pe and sr
# - Most variants (~80%) are complex rearrangements. Limiting coverage does not affect the fraction of variant of each class, maybe with the exception of duplications (DUP): there are 7 times more duplications when no regions are excluded based on their coverage, compared to the ~3-fold average across all types --this probably reflects the fact that regions with high coverage reflect a good number of true duplications

# In[295]:

df = svs['no_coverage_limit']


# In[296]:

df[df['TYPE'] == 'BND'].head()


# In[399]:

colors()[0]


# In[400]:

28821/8277.


# In[ ]:



