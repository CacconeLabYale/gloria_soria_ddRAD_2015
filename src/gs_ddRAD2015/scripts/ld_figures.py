
# coding: utf-8

# # Purpose:
# 
# 2015-03-13 (Friday)
# 
# Explore and characterize the results of the learned Beta filter method.

# # Contents
# [Loading-files](#Loading-files)
# 
# [Filtering-Stats](#Filtering-Stats)
# 
# [Characterization-of-contigs-with/without-regard-to-selected-SNP-pairs](#Characterization-of-contigs-with/without-regard-to-selected-SNP-pairs)
# - [All-contigs](#All-contigs)
# - [All-SNP-pair-contigs](#All-SNP-pair-contigs)
# - [LD-filtered-SNP-pair-contigs](#LD-filtered-SNP-pair-contigs)
# 
# [Compare-with-Tajima's-D-results-from-AndreaG](#Compare-with-Tajima's-D-results-from-AndreaG)
# - [Examine-ld_contig_taj_win_filter](#Examine-ld_contig_taj_win_filter)
# - [Basic-summary-table-referencing-SNP-pairs](#Basic-summary-table-referencing-SNP-pairs)
# - [How-is-q-value-related-to-SNP-pair-distance?](#How-is-q-value-related-to-SNP-pair-distance?)
# - [What-is-the-distribution-of-SNPs-per-bin-used-for-Tajima's-D-bin_50](#What-is-the-distribution-of-SNPs-per-bin-used-for-Tajima's-D-bin_50)
# - [How-Tajima's-D-score-bin_50-relate-to-number-of-SNPs-in-the-bin](#How-Tajima's-D-score-bin_50-relate-to-number-of-SNPs-in-the-bin)
# 
# [Average-LD-per-bin](#Average-LD-per-bin)
# - [All-contigs-together](#All-contigs-together)
# - [5-random-contigs](#5-random-contigs)
# - 

# ## Imports:

# In[3]:

import datetime as dt

import matplotlib.pyplot as plt
import seaborn as sns
import ggplot as gp


import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 60)
# import tables as h5

import itertools as it
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy
from scikits import bootstrap as bs
import statsmodels.api as sm
import statsmodels.stats.multitest as smm

import munch

import pymc as mc

from spartan.utils.genome_specific.GfusI1 import GfusI1_0
from spartan.utils.fastas import ParseFastA


# In[4]:

# set figure characteristics

# size
sns.set_context("talk")

# style
sns.set_style("whitegrid")

# ## File paths:

# In[5]:

# define paths to files

base_out_dir = "/home/gus/Documents/YalePostDoc/project_stuff/g_f_fucipes_uganda/ddrad58/manuscript/figures/ld"

contig_name_length_path = "/home/gus/Dropbox/uganda_data/data_repos/genome_info/assembly_info/contig_name_length.csv"

ld_results_pickle="/home/gus/Documents/YalePostDoc/project_stuff/g_f_fucipes_uganda/ddrad58/ld_thresholds/post_MAP_calc.plk"
tajimas_csv = "/home/gus/Documents/YalePostDoc/project_stuff/g_f_fucipes_uganda/ddrad58/data_from_andrea/Tajima50.csv"



# # Helper functions

###########################

def recode_taj_chrom(df):
    recode_func = lambda x: x.split(':')[-1]

    CHROM = df.CHROM.apply(recode_func)
    df.CHROM = CHROM


# # Loading files

# In[10]:

# load our results tables
ld = pd.read_pickle(ld_results_pickle)
ld.head()


# In[11]:

contig_info = pd.read_csv(contig_name_length_path)
contig_info.head()


# In[12]:

taj50 = pd.read_csv(tajimas_csv, sep='\t')
taj50.head()


# In[13]:

recode_taj_chrom(taj50)
taj50.head()


# # Filtering Stats

# ## SNP-pairs in all bins at BH corrected $p \leq 0.01$

# In[14]:

sum(ld.one_minus_cdf_BH <= 0.01)


# ## SNP-pairs in all bins at BH corrected $p \le 0.05$

# In[15]:

sum(ld.one_minus_cdf_BH <= 0.05)


# ## Lowest $r^2$ retained at  $p \le 0.05$ or $0.01$

# In[16]:

q_05 = ld.query("one_minus_cdf_BH <= 0.05")
q_05.R2.min()


# In[17]:

q_01 = ld.query("one_minus_cdf_BH <= 0.01")
q_01.R2.min()


# ## How many SNP-pairs have  $r^2 \ge 0.82$?

# In[18]:

sum(ld.R2 >= 0.82)


# In[19]:

1-(5284.0/26495)


# ## Characterization of contigs with/without regard to selected SNP-pairs

# In[20]:

# join contig length and kk_name contig info to the LD table
ld_contig = pd.merge(left=ld, right=contig_info, how='inner', left_on="CHR_A", right_on="scaf_name")
ld_contig.head()


# ### All contigs

# In[21]:

sns.distplot(contig_info.length, color="coral", kde=0);
median_contig_len = contig_info.length.median()
mean_contig_len = contig_info.length.mean()
plt.text(median_contig_len,400,"*median length = {}".format(median_contig_len), fontsize=14);
plt.text(mean_contig_len,200,"*mean length = {}".format(mean_contig_len), fontsize=14);
plt.text(mean_contig_len,800,"Number of contigs = {}".format(len(contig_info)), fontsize=14);
plt.title("All contigs");


# ### All SNP-pair contigs

# In[22]:

sp_contigs = contig_info[contig_info.scaf_name.isin(ld_contig.scaf_name.unique())]
len(sp_contigs)


# In[23]:

sns.distplot(sp_contigs.length, color="coral", kde=0);
sp_median_contig_len = sp_contigs.length.median()
sp_mean_contig_len = sp_contigs.length.mean()
plt.text(sp_median_contig_len,100,"*median length = {}".format(sp_median_contig_len), fontsize=14);
plt.text(sp_mean_contig_len,50,"*mean length = {}".format(sp_mean_contig_len), fontsize=14);
plt.text(sp_mean_contig_len,200,"Number of contigs = {}".format(len(sp_contigs)), fontsize=14);
plt.title("Contigs with a SNP-pair");


# ### LD filtered SNP-pair contigs

# In[24]:

ld_filt_contigs = ld_contig.query("one_minus_cdf_BH <= 0.01")
ld_filt_contigs = contig_info[contig_info.scaf_name.isin(ld_filt_contigs.scaf_name.unique())]
len(ld_filt_contigs)


# In[25]:

sns.distplot(ld_filt_contigs.length, color="coral", kde=0);
ld_filt_median_contig_len = ld_filt_contigs.length.median()
ld_filt_mean_contig_len = ld_filt_contigs.length.mean()
plt.text(ld_filt_median_contig_len, 40, "*median length = {}".format(ld_filt_median_contig_len), fontsize=14);
plt.text(ld_filt_mean_contig_len, 20, "*mean length = {}".format(ld_filt_mean_contig_len), fontsize=14);
plt.text(ld_filt_mean_contig_len ,80, "Number of contigs = {}".format(len(ld_filt_contigs)), fontsize=14);
plt.title(r"Contigs with a SNP-pair filtered by binned LD ($q \leq 0.01$)");


# ## Compare with Tajima's D results from AndreaG

# In[26]:

taj50.head()


# In[27]:

ld_contig.head()


# ### To accomplish
# Need to filter out `ld_contig` data that has either SNP ocurring in the bins defined by `taj50.CHROM:taj50.BIN_start-[taj50.BIN_start+50]`
# 
# First try:
# 
# - join INNER `ld_contig` and `taj50` on `left_on=kk_name`, `right_on=CHROM` as `ld_contig_taj`
# - reatain those rows where `ld_contig_taj.BP_A` or `ld_contig_taj.BP_B` is inside `ld_contig_taj.BIN_start` - `ld_contig_taj.BIN_start+50`

# In[28]:

ld_contig_taj = pd.merge(left=ld_contig, right=taj50, how='inner', left_on='kk_name', right_on='CHROM')


# In[29]:

def get_taj_win_mask(df):
    taj_win_start = df.BIN_start
    taj_win_end = df.BIN_start + 50
    
    a_mask = (ld_contig_taj.BIN_start <= ld_contig_taj.BP_A) & (ld_contig_taj.BP_A <= ld_contig_taj.BIN_start + 50)
    b_mask = (ld_contig_taj.BIN_start <= ld_contig_taj.BP_B) & (ld_contig_taj.BP_B <= ld_contig_taj.BIN_start + 50)
    
    return (a_mask | b_mask)

# get all at first
ld_contig_taj_win = ld_contig_taj[get_taj_win_mask(ld_contig_taj)]

# now subset these to only rows that meet the LD bin filter
ld_contig_taj_win_filter = ld_contig_taj_win.query("one_minus_cdf_BH <= 0.01")


# In[30]:

print len(ld_contig_taj)
print len(ld_contig_taj_win)
print len(ld_contig_taj_win_filter)


# In[31]:

ld_contig_taj_win_filter.head(15)


# ## Examine `ld_contig_taj_win_filter` 

# ### Basic summary table referencing SNP-pairs

# In[32]:

ld_contig_taj_win_filter_t1 = pd.pivot_table(ld_contig_taj_win_filter,
                                             index=['scaf_name','BP_A','BP_B','distance_bin'],
                                             fill_value=0,
                                            )
ld_contig_taj_win_filter_t1.head()


# ### How is q-value related to SNP-pair distance?

# In[33]:

sns.jointplot(x='BP_DELTA',
              y='one_minus_cdf_BH', 
              data=ld_contig_taj_win_filter, 
              kind='reg',
              color='lightblue',
              xlim=(0,ld_contig_taj_win_filter.distance_bin.max()),
              ylim=(0,ld_contig_taj_win_filter.one_minus_cdf_BH.max()));


# ### What is the distribution of SNPs per bin used for Tajima's D bin_50

# In[34]:

sns.distplot(ld_contig_taj_win_filter.N_SNPs, color='lightblue', kde=0);


# ### How Tajima's D score bin_50 relate to number of SNPs in the bin

# In[35]:

sns.jointplot(x='N_SNPs',
              y='TajimaD', 
              data=ld_contig_taj_win_filter, 
              kind='reg',
              color='lightblue',
             );


# # Average LD per bin

# ### How many contigs are available to each distance_bin?

# In[36]:

def get_contigs_per_bin(d_bins,contig_info):
    
    cpb = {}
    
    for b in d_bins:
        cpb[b] = sum(contig_info.length > b)
        
    return pd.Series(cpb)


# In[37]:

d_bins = ld_contig.distance_bin.unique()
d_bins.sort()
d_bins
# Generate a dict to map how many contigs avail to each bin
contigs_per_bin = get_contigs_per_bin(d_bins,contig_info)
contigs_per_bin


# In[38]:

contigs_per_bin = pd.DataFrame(contigs_per_bin,columns=['contigs_per_bin'])


# In[39]:

contigs_per_bin = contigs_per_bin.reset_index().rename(columns={'index':'distance_bin'}, inplace=False)


# In[40]:

contigs_per_bin.head()


# In[41]:

gp.ggplot(gp.aes(x='distance_bin', y='contigs_per_bin'), data=contigs_per_bin) +     gp.geom_point(color='lightblue') +     gp.stat_smooth(span=.15, color='red', se=True)


# ## All contigs together

# In[42]:

ld_contig.head()


# In[43]:

mean_bin_r2_all = ld_contig.groupby("distance_bin").mean().reset_index()


# In[44]:

median_bin_r2_all = ld_contig.groupby("distance_bin").median().reset_index()


# In[45]:

gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=mean_bin_r2_all.query("distance_bin >= 0")) +     gp.geom_point(color='lightblue', alpha=0.4) +     gp.stat_smooth(span=.15, color='red', se=True)


# In[46]:

gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=mean_bin_r2_all.query("distance_bin <= 30000")) +     gp.geom_point(color='lightblue') +     gp.stat_smooth(span=0.15, color='red', se=True) 


# In[47]:

len(ld_contig)


# In[48]:

len_contigs_per_bin = ld_contig.pivot_table(index=['distance_bin'], 
                                            values=['scaf_name'],
                                            aggfunc=[len]
                                           )['len'].reset_index()

len_contigs_per_bin = len_contigs_per_bin.rename(columns={'scaf_name':'SNP-pairs'}, inplace=False)
len_contigs_per_bin.head()


# In[49]:

len_contigs_per_bin.


# In[58]:

sp_gt100mask = len_contigs_per_bin["SNP-pairs"] > 100

gp.ggplot(gp.aes(x='distance_bin', y='SNP-pairs'), data=len_contigs_per_bin[sp_gt100mask]) +     gp.geom_point(color='lightblue') +     gp.stat_smooth(span=.15, color='red', se=True) + gp.ylab('SNP-pairs per bin')


# In[57]:




# In[56]:

ld_contig.head()


# In[57]:

d_bin_v_others = ld_contig.pivot_table(index=['distance_bin'], 
                                        values=['R2','one_minus_cdf_BH'],
                                        aggfunc=[np.mean]
                                       )['mean'].reset_index()
d_bin_v_others.head()


# In[58]:

d_bin_v_others = d_bin_v_others.merge(right=len_contigs_per_bin, 
                     how='inner', 
                     on='distance_bin'
                     ).merge(right=contigs_per_bin, 
                             how='inner', 
                             on='distance_bin'
                            )
d_bin_v_others.head()


# In[59]:

d_bin_v_others_melt = pd.melt(d_bin_v_others, id_vars=['distance_bin'])


# In[81]:

d_bin_v_others_melt.head()


# In[84]:

xmin = 0
xmax = 100000000000
gp.ggplot(gp.aes(x='distance_bin', y='value'), 
          data=d_bin_v_others_melt.query("{xmin} <= distance_bin <= {xmax}".format(xmin=xmin,
                                                                                   xmax=xmax
                                                                                  ))) + \
    gp.geom_point(color='lightblue', alpha=0.06) + \
    gp.stat_smooth(span=0.2, color='red', se=True)  + \
    gp.facet_wrap("variable")


# In[61]:

# g = sns.PairGrid(d_bin_v_others.loc[:,["R2", "one_minus_cdf_BH",    "SNP-pairs",   "contigs_per_bin",]])
# g.map_upper(plt.scatter)
# g.map_lower(sns.kdeplot, cmap="Blues_d")
# g.map_diag(sns.kdeplot, lw=3, legend=False)
d_bin_vars = d_bin_v_others.loc[:,["R2", "one_minus_cdf_BH",    "SNP-pairs",   "contigs_per_bin",]]


# In[62]:

d_bin_vars.head()


# In[63]:

sns.palplot(sns.cubehelix_palette(8, start=.5, rot=-.75))


# In[64]:

my_cmap=sns.cubehelix_palette(40, start=.5, rot=-.75, as_cmap=True)
cc = sns.mpl.colors.ColorConverter()
marginal_color = cc.to_rgb(arg=my_cmap.colors[int(255*1)])


# In[65]:

# sns.jointplot('SNP-pairs','R2',d_bin_vars, kind='kde',
#               joint_kws=dict(shade=True,
#                              cmap=my_cmap,
#                              n_levels=40
#                             ),
#               marginal_kws=dict(shade=True, color=my_cmap.colors[int(256*0.66)])
#              )

g = sns.JointGrid('SNP-pairs','R2',d_bin_vars)
g.plot_marginals(sns.distplot, kde=False, color=marginal_color)
g.plot_joint(sns.kdeplot, shade=True, cmap=my_cmap, n_levels=40);


# In[66]:

# sns.jointplot('one_minus_cdf_BH','R2',d_bin_vars, kind='kde',
#               joint_kws=dict(shade=True,
#                              cmap=my_cmap,
#                              n_levels=40
#                             ),
#               marginal_kws=dict(shade=True, color=my_cmap.colors[int(256*0.66)])
#              )

g = sns.JointGrid('one_minus_cdf_BH','R2',d_bin_vars)
g.plot_marginals(sns.distplot, kde=False, color=marginal_color)
g.plot_joint(sns.kdeplot, shade=True, cmap=my_cmap, n_levels=40, alpha=1);


# In[70]:

sns.jointplot('contigs_per_bin','SNP-pairs',d_bin_vars, kind='kde',
              joint_kws=dict(shade=True,
                             cmap=my_cmap,
                             n_levels=8
                            ),
              marginal_kws=dict(shade=True, color=my_cmap.colors[int(256*0.66)])
             )


# ## 5 random contigs
