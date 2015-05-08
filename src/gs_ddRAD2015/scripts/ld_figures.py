import ipdb



import inspect

import collections
import os
import matplotlib.pyplot as plt
import munch
import seaborn as sns
import ggplot as gp

import pandas as pd

pd.set_option('display.max_columns', 60)

import numpy as np


def run(ld_pickle, out_dir, formats, contig_length):
    """
    Takes processed LD d in form of a python pickle and produces figures.

    Args:
        ld_pickle (str): path to pickle file.
        out_dir (str): path to directory where the figures should go.
        formats (list): one or more of ['png','svg','pdf','none'].
        contig_length (str): path to csv file with two labeled columns: ['scaf_name','length'].

    Returns:
        None
    """

    figs = Figures(out_dir, formats, extras=False)



    # style
    sns.set_context("talk")
    sns.set_style("whitegrid")

    # # Loading files
    # load our results tables
    figs.d.ld = ld = pd.read_pickle(ld_pickle)

    figs.d.contig_info = contig_info = pd.read_csv(contig_length)

    # join contig length and contig info to the LD table
    figs.d.ld_contig = ld_contig = pd.merge(left=ld, right=contig_info, how='inner', left_on="CHR_A", right_on="scaf_name")

    # list all contigs that had at least one Snp-pair
    figs.d.sp_contigs = sp_contigs = contig_info[contig_info.scaf_name.isin(ld_contig.scaf_name.unique())]

    ipdb.set_trace()
    # get list of distance_bins
    d_bins = ld_contig.distance_bin.unique()
    d_bins.sort()
    figs.d.d_bins = d_bins

    # Generate a dict to map how many contigs avail to each bin
    contigs_per_bin = self.get_contigs_per_bin(d_bins, contig_info)
    contigs_per_bin = pd.DataFrame(contigs_per_bin, columns=['contigs_per_bin'])
    figs.d.contigs_per_bin = contigs_per_bin = contigs_per_bin.reset_index().rename(columns={'index': 'distance_bin'}, inplace=False)

    figs.d.mean_bin_r2_all = mean_bin_r2_all = ld_contig.groupby("distance_bin").mean().reset_index()
    # median_bin_r2_all = ld_contig.groupby("distance_bin").median().reset_index()

    ##############################################
    ipdb.set_trace()

    len_contigs_per_bin = ld_contig.pivot_table(index=['distance_bin'],
                                                values=['scaf_name'],
                                                aggfunc=[len]
                                                )['len'].reset_index()

    figs.d.len_contigs_per_bin = len_contigs_per_bin = len_contigs_per_bin.rename(columns={'scaf_name': 'SNP-pairs'}, inplace=False)

    ###############################################

    d_bin_v_others = ld_contig.pivot_table(index=['distance_bin'],
                                           values=['R2', 'one_minus_cdf_BH'],
                                           aggfunc=[np.mean]
                                           )['mean'].reset_index()

    figs.d.d_bin_v_others = d_bin_v_others = d_bin_v_others.merge(right=len_contigs_per_bin,
                                          how='inner',
                                          on='distance_bin'
                                          ).merge(right=contigs_per_bin,
                                                  how='inner',
                                                  on='distance_bin'
                                                  )

    figs.d.d_bin_v_others_melt = d_bin_v_others_melt = pd.melt(d_bin_v_others, id_vars=['distance_bin'])

    #################################################

    figs.self.make_figures()

# ######################################################################



# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################




class Figures(object):

    def __init__(self, out_dir, formats, save_plots=False, extras=False):
        self.save_plots = save_plots
        self.extras = extras
        self.base_dir = out_dirgitzp
        self.formats = formats
        self.d = munch.Munch()
        self.runners = self._get_runners()

    def save_figs(self, base_dir, fname, save_types, is_ggplot=False):
        assert isinstance(base_dir, str)
        assert isinstance(fname, str)
        assert isinstance(save_types, collections.Iterable)
        assert (is_ggplot is False) or isinstance(is_ggplot, gp.ggplot)

        if not self.save_plots:
            print "WARNING: 'save_plots' set to False."
            return None

        path = "{pth}.{{ext}}".format(pth=os.path.join(base_dir, fname))

        for t in save_types:

            if is_ggplot:
                gp.ggsave(path.format(ext=t), plot=is_ggplot)
            else:
                plt.savefig(path.format(ext=t), bbox_inches='tight')
                print "Saved {0}.".format(path.format(ext=t))

    def plot_bin_dists(self, df, bin_def="distance_bin <= 500"):
        g = sns.FacetGrid(data=df.query(bin_def),
                          col="distance_bin",
                          sharey=False,
                          sharex=False,
                          col_wrap=4,
                          )

        return g.map(plt.hist, "R2", color="coral")

    def get_contigs_per_bin(self, d_bins, contig_info):
        cpb = {}

        for b in d_bins:
            cpb[b] = sum(contig_info.length > b)

        return pd.Series(cpb)

    def plot_scat_w_line(self, gp_aes):
        return gp_aes + gp.geom_point(color='coral') + gp.stat_smooth(span=.2, color='blue',
                                                                      se=False) + gp.theme_seaborn(
            context='talk')

    def make_figures(self):

        for name, runner in self.runners:
            if not self.extras:
                if name != 'display_extras':
                    runner()
            else:
                runner()

    def _get_runners(self):
        runners = inspect.getmembers(self, lambda x: inspect.ismethod(x))

        return runners

    # ########## Figure  ##########
    def distance_bins_0_500(self):
        fname = "distance_bins_0_500"

        p = self.plot_bin_dists(ld, bin_def="distance_bin <= 500")

        p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )

    # ########## Figure  ##########
    def distance_bins_7000_7500(self):
        fname = "distance_bins_7000_7500"

        p = self.plot_bin_dists(ld, bin_def="7000 <= distance_bin <= 7500")

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )

    # ########## Figure  ##########
    def distance_bins_15000_15500(self):
        fname = "distance_bins_15000_15500"

        p = self.plot_bin_dists(ld, bin_def="15000 <= distance_bin <= 15500")

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )

    # ########## Figure  ##########
    def distance_bins_30000_30500(self):
        fname = "distance_bins_30000_30500"

        p = self.plot_bin_dists(ld, bin_def="30000 <= distance_bin <= 30500")

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )

    # ########## Figure  ##########
    def distance_bins_50000_50500(self):
        fname = "distance_bins_50000_50500"

        p = self.plot_bin_dists(ld, bin_def="50000 <= distance_bin <= 50500")

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )

    # ########## Figure  ##########
    def distance_bins_70000_70500(self):
        fname = "distance_bins_70000_70500"

        p = self.plot_bin_dists(ld, bin_def="70000 <= distance_bin <= 70500")

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )

    # ########## Figure  ##########
    def distance_bins_100000_100500(self):
        fname = "distance_bins_100000_100500"

        p = self.plot_bin_dists(ld, bin_def="100000 <= distance_bin <= 100500")

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )

    # ########## Figure  ##########
    def all_contig_len_dist(self):
        fname = "all_contig_len_dist"

        sns.distplot(contig_info.length, color="coral", kde=0);
        median_contig_len = contig_info.length.median()
        mean_contig_len = contig_info.length.mean()
        plt.text(median_contig_len, 400, "*median length = {}".format(median_contig_len), fontsize=14);
        plt.text(mean_contig_len, 200, "*mean length = {}".format(mean_contig_len), fontsize=14);
        plt.text(mean_contig_len, 800, "Number of contigs = {}".format(len(contig_info)), fontsize=14);
        plt.title("All contigs");

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=('png', 'pdf', 'svg')
                  )

    # ### All SNP-pair contigs
    # ########## Figure  ##########
    def all_sp_contig_len_dist(self):
        ipdb.set_trace()
        fname = "all_sp_contig_len_dist"

        sns.distplot(sp_contigs.length, color="coral", kde=0);
        sp_median_contig_len = sp_contigs.length.median()
        sp_mean_contig_len = sp_contigs.length.mean()
        plt.text(sp_median_contig_len, 100, "*median length = {}".format(sp_median_contig_len), fontsize=14);
        plt.text(sp_mean_contig_len, 50, "*mean length = {}".format(sp_mean_contig_len), fontsize=14);
        plt.text(sp_mean_contig_len, 200, "Number of contigs = {}".format(len(sp_contigs)), fontsize=14);
        plt.title("Contigs with a SNP-pair");

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=('png', 'pdf', 'svg')
                  )

    # ########## Figure  ##########
    def distance_VS_contigs_per_bin(self):
        fname = "distance_VS_contigs_per_bin"

        p = self.plot_scat_w_line(gp.ggplot(gp.aes(x='distance_bin', y='contigs_per_bin'), data=contigs_per_bin))

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )

    # ## All contigs together


    # ########## Figure  ##########
    def distance_VS_r2_all(self):
        ipdb.set_trace()
        fname = "distance_VS_r2_all"

        p = self.plot_scat_w_line(
            gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=mean_bin_r2_all.query("distance_bin >= 0")))

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )

    # ########## Figure  ##########
    def distance_VS_r2_le30K(self):
        fname = "distance_VS_r2_le30K"

        p = self.plot_scat_w_line(
            gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=mean_bin_r2_all.query("distance_bin <= 30000")))

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )

    # ########## Figure  ##########
    def distance_VS_r2_le40K(self):
        fname = "distance_VS_r2_le40K"

        p = self.plot_scat_w_line(
            gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=mean_bin_r2_all.query("distance_bin <= 40000")))

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )

    # ########## Figure  ##########
    def distance_VS_r2_le50K(self):
        fname = "distance_VS_r2_le50K"

        p = self.plot_scat_w_line(
            gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=mean_bin_r2_all.query("distance_bin <= 50000")))

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )

    # ########## Figure  ##########
    def distance_VS_sp_per_bin_gt100(self):
        fname = "distance_VS_sp_per_bin_gt100"

        sp_gt100mask = len_contigs_per_bin["SNP-pairs"] > 100

        p = self.plot_scat_w_line(gp.ggplot(gp.aes(x='distance_bin',
                                              y='SNP-pairs'),
                                       data=len_contigs_per_bin[sp_gt100mask])
                             ) + gp.ylab('SNP-pairs per bin')

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )

    # ########## Figure  ##########
    def distance_VS_sp_per_bin_gt150(self):
        fname = "distance_VS_sp_per_bin_gt150"

        sp_gt150mask = len_contigs_per_bin["SNP-pairs"] > 150

        p = self.plot_scat_w_line(gp.ggplot(gp.aes(x='distance_bin',
                                              y='SNP-pairs'),
                                       data=len_contigs_per_bin[sp_gt150mask])
                             ) + gp.ylab('SNP-pairs per bin')

        print p

        self.save_figs(base_dir=self.d.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )

    # In[ ]:
    def distance_VS_avgR2_snpperbin_contigsperbin_q_b_bmin_to_bbmax(self):
        extents = ((0, 1000),
                   (0, 3000),
                   (150, 1000),
                   (150, 3000),
                   (150, 10000),
                   (150, 20000),
                   (150, 30000),
                   (150, 40000),
                   (150, 50000),
                   (150, 60000),
                   (150, 70000),
                   (150, 80000),
                   (150, 90000),
                   (150, 100000),)

        xmins = (0, 150,)
        xmaxs = (1000,)

        for xmin, xmax in extents:
            fname = "distance_VS_avgR2_snpperbin_contigsperbin_q_b{bmin}-to-b{bmax}".format(
                bmin=xmin,
                bmax=xmax
            )

            p = self.plot_scat_w_line(
                gp.ggplot(
                    gp.aes(
                        x='distance_bin', y='value'),
                    data=d_bin_v_others_melt.query(
                        "{xmin} <= distance_bin <= {xmax}".format(
                            xmin=xmin,
                            xmax=xmax
                        )
                    )
                )
            ) + \
                gp.facet_wrap("variable") + \
                gp.ggtitle(fname)

            print p

            self.save_figs(base_dir=self.d.base_dir,
                      fname=fname,
                      save_types=self.formats,
                      is_ggplot=p
                      )

    def display_extras(self):
        # In[ ]:

        # g = sns.PairGrid(d_bin_v_others.loc[:,["R2", "one_minus_cdf_BH",    "SNP-pairs",   "contigs_per_bin",]])
        # g.map_upper(plt.scatter)
        # g.map_lower(sns.kdeplot, cmap="Blues_d")
        # g.map_diag(sns.kdeplot, lw=3, legend=False)
        d_bin_vars = d_bin_v_others.loc[:, ["R2", "one_minus_cdf_BH", "SNP-pairs", "contigs_per_bin", ]]

        sns.palplot(sns.cubehelix_palette(8, start=.5, rot=-.75))


        # In[ ]:

        my_cmap = sns.cubehelix_palette(40, start=.5, rot=-.75, as_cmap=True)
        cc = sns.mpl.colors.ColorConverter()
        marginal_color = cc.to_rgb(arg=my_cmap.colors[int(255 * 1)])


        # In[ ]:

        # sns.jointplot('SNP-pairs','R2',d_bin_vars, kind='kde',
        #               joint_kws=dict(shade=True,
        #                              cmap=my_cmap,
        #                              n_levels=40
        #                             ),
        #               marginal_kws=dict(shade=True, color=my_cmap.colors[int(256*0.66)])
        #              )

        g = sns.JointGrid('SNP-pairs', 'R2', d_bin_vars)
        g.plot_marginals(sns.distplot, kde=False, color=marginal_color)
        g.plot_joint(sns.kdeplot, shade=True, cmap=my_cmap, n_levels=40);


        # In[ ]:

        # sns.jointplot('one_minus_cdf_BH','R2',d_bin_vars, kind='kde',
        #               joint_kws=dict(shade=True,
        #                              cmap=my_cmap,
        #                              n_levels=40
        #                             ),
        #               marginal_kws=dict(shade=True, color=my_cmap.colors[int(256*0.66)])
        #              )

        g = sns.JointGrid('one_minus_cdf_BH', 'R2', d_bin_vars)
        g.plot_marginals(sns.distplot, kde=False, color=marginal_color)
        g.plot_joint(sns.kdeplot, shade=True, cmap=my_cmap, n_levels=40, alpha=1);


        # In[ ]:

        sns.jointplot('contigs_per_bin', 'SNP-pairs', d_bin_vars, kind='kde',
                      joint_kws=dict(shade=True,
                                     cmap=my_cmap,
                                     n_levels=8
                                     ),
                      marginal_kws=dict(shade=True, color=my_cmap.colors[int(256 * 0.66)])
                      )





