import click
import ipdb
import sh



import inspect
import collections
import os
import matplotlib
# turn off interactive plotting.  It kills memory to keep all these figs open.
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import munch
import seaborn as sns
import ggplot as gp
import pandas as pd

pd.set_option('display.max_columns', 60)

import numpy as np

echo = click.echo

# def echo(msg):
#     click.echo(msg)
#     sh.fembot(t=msg)





def run(ld_pickle, out_dir, formats, contig_length, save_tables, force_save):
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

    skip_calcs = False
    data_dir = os.path.join(out_dir, 'data_tables')
    if save_tables:
        if not os.path.exists(data_dir):
            echo("\nld_figures: No data_tables dir found, NOT skipping data crunching stage.")
            figs = Figures(out_dir, formats, extras=False)
        else:
            try:
                echo("\nld_figures: Attempting to load data_tables.")
                figs = Figures(out_dir, formats, extras=False)
                figs.load_data_tables()
                skip_calcs = True
            except IOError as exc:
                echo("\nld_figures: failed to load one or more files, NOT skipping data crunching stage.")
                print(exc)
                skip_calcs = False

    else:
        figs = Figures(out_dir, formats, extras=False)

    # style
    sns.set_context("talk")
    sns.set_style("whitegrid")

    if skip_calcs:
        echo("\nld_figures: Skipping data crunching stage due to found data_tables.")

    else:
        echo("\nld_figures: Entering data crunching stage:\nld_figures: This may take a while..."
              # "but please enjoy this fancy progress bar!"
              "\nld_figures: Really. Maybe you should go get some coffee...")

        # with click.progressbar(length=12) as progbar:

        # # Loading files
        # load our results tables
        figs.d.ld = pd.read_pickle(ld_pickle)
        # progbar.update(1)

        figs.d.contig_info = pd.read_csv(contig_length)
        # progbar.update(2)

        # join contig length and contig info to the LD table
        figs.d.ld_contig = pd.merge(left=figs.d.ld, right=figs.d.contig_info, how='inner', left_on="CHR_A",
                                    right_on="scaf_name")
        # progbar.update(3)

        # list all contigs that had at least one Snp-pair
        figs.d.sp_contigs = figs.d.contig_info[figs.d.contig_info.scaf_name.isin(figs.d.ld_contig.scaf_name.unique())]
        # progbar.update(4)
        #
        # ipdb.set_trace()
        # get list of distance_bins
        d_bins = figs.d.ld_contig.distance_bin.unique()
        d_bins.sort()
        figs.d.d_bins = pd.Series(d_bins)
        # progbar.update(5)
        #
        # Generate a dict to map how many contigs avail to each bin
        contigs_per_bin = figs.get_contigs_per_bin(d_bins, figs.d.contig_info)
        contigs_per_bin = pd.DataFrame(contigs_per_bin, columns=['contigs_per_bin'])
        figs.d.contigs_per_bin = contigs_per_bin.reset_index().rename(columns={'index': 'distance_bin'}, inplace=False)
        # progbar.update(6)
        #
        figs.d.mean_bin_r2_all = figs.d.ld_contig.groupby("distance_bin").mean().reset_index()
        # median_bin_r2_all = ld_contig.groupby("distance_bin").median().reset_index()
        # progbar.update(7)
        #
        ##############################################
        # ipdb.set_trace()

        figs.d.len_contigs_per_bin = figs.d.ld_contig.pivot_table(index=['distance_bin'],
                                                    values=['scaf_name'],
                                                    aggfunc=[len]
                                                    )['len'].reset_index()
        # progbar.update(8)

        figs.d.len_contigs_per_bin = figs.d.len_contigs_per_bin.rename(columns={'scaf_name': 'SNP-pairs'}, inplace=False)
        # progbar.update(9)
        ###############################################

        figs.d.d_bin_v_others = figs.d.ld_contig.pivot_table(index=['distance_bin'],
                                               values=['R2', 'one_minus_cdf_BH'],
                                               aggfunc=[np.mean]
                                               )['mean'].reset_index()
        # progbar.update(10)

        figs.d.d_bin_v_others = figs.d.d_bin_v_others.merge(right=figs.d.len_contigs_per_bin,
                                                                      how='inner',
                                                                      on='distance_bin'
                                                                      ).merge(right=figs.d.contigs_per_bin,
                                                                              how='inner',
                                                                              on='distance_bin'
                                                                              )
        # progbar.update(11)

        figs.d.d_bin_v_others_melt = pd.melt(figs.d.d_bin_v_others, id_vars=['distance_bin'])
        # progbar.update(12)

        if save_tables and (force_save or not skip_calcs):
            figs.write_data_tables()
        else:
            pass

    #################################################
    echo("\nld_figures: Data crunching done. Beginning plotting stage.")
    figs.make_figures()

# ######################################################################



# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################




class Figures(object):

    def __init__(self, out_dir, formats, extras=False):

        self.extras = extras
        self.base_dir = out_dir
        self.data_tables = os.path.join(out_dir, 'data_tables')
        self.formats = formats
        self.d = munch.Munch()
        self.runners = None
        self.force_write = False
        self.d_keys = ['ld',
                       'contig_info',
                       'contigs_per_bin',
                       'd_bin_v_others',
                       'd_bin_v_others_melt',
                       'd_bins',
                       'ld_contig',
                       'len_contigs_per_bin',
                       'mean_bin_r2_all',
                       'sp_contigs']

    def load_data_tables(self):
        template = "{dir}/{key}"

        echo("\nld_figures: loading data_tables.")
        with click.progressbar(self.d_keys) as keys:
            for key in keys:
                self.d[key] = pd.read_csv(template.format(dir=self.data_tables, key=key),
                                          sep=',',
                                          na_values=['no_val'],
                                          keep_default_na=False)

    def write_data_tables(self):
        template = "{dir}/{key}"

        echo("\nld_figures: writing data_tables.")
        with click.progressbar(self.d.items()) as items:
            for key, value in items:
                opath = template.format(dir=self.data_tables, key=key)
                if not os.path.exists(opath) or self.force_write:
                    try:
                        value.to_csv(opath,
                                     sep=',',
                                     na_rep='no_val')

                    except IOError as exc:
                        if "No such file or directory" in exc.args[-1]:
                            os.mkdir(self.data_tables)

                            value.to_csv(opath,
                                         sep=',',
                                         na_rep='no_val')

    def save_figs(self, base_dir, fname, save_types, is_ggplot=False):
        assert isinstance(base_dir, str)
        assert isinstance(fname, str)
        assert isinstance(save_types, collections.Iterable)
        assert (is_ggplot is False) or isinstance(is_ggplot, gp.ggplot)

        if 'none' in self.formats:
            echo("\nld_figures: WARNING: 'none' in --formats, figures are NOT being saved.")
            return None

        path = "{pth}.{{ext}}".format(pth=os.path.join(base_dir, fname))

        for t in save_types:
            try:
                if is_ggplot:
                    gp.ggsave(path.format(ext=t), plot=is_ggplot)
                    click.echo("\nld_figures: Saved {0}.".format(path.format(ext=t)))
                else:
                    plt.savefig(path.format(ext=t), bbox_inches='tight')
                    click.echo("\nld_figures: Saved {0}.".format(path.format(ext=t)))
            except IndexError as exc:
                if 'index out of bounds' in exc.args[0]:
                    click.echo("\nld_figures: skipping due to lack of data.")

    def _plot_bin_dists(self, df, bin_def="distance_bin <= 500"):
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

    def _plot_scat_w_line(self, gp_aes):
        return gp_aes + gp.geom_point(color='coral') + gp.stat_smooth(span=.2, color='blue',
                                                                      se=False) + gp.theme_seaborn(
            context='talk')

    def make_figures(self):
        self._set_runners()

        with click.progressbar(self.runners) as runners:
            for name, runner in runners:

                if name.startswith('plot_'):
                    if not self.extras:
                        if name != 'display_extras':
                            runner()
                            # plt.show()
                    else:
                        runner()
                        # plt.show()

    def _set_runners(self):
        runners = inspect.getmembers(self, lambda x: inspect.ismethod(x))

        self.runners = runners

    # ########## Figure  ##########
    def plot_distance_bins_0_500(self):
        fname = "distance_bins_0_500"

        p = self._plot_bin_dists(self.d.ld, bin_def="distance_bin <= 500")

        p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )

    # ########## Figure  ##########
    def plot_distance_bins_7000_7500(self):
        fname = "distance_bins_7000_7500"

        p = self._plot_bin_dists(self.d.ld, bin_def="7000 <= distance_bin <= 7500")

        # # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )
        

    # ########## Figure  ##########
    def plot_distance_bins_15000_15500(self):
        fname = "distance_bins_15000_15500"

        p = self._plot_bin_dists(self.d.ld, bin_def="15000 <= distance_bin <= 15500")

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )
        

    # ########## Figure  ##########
    def plot_distance_bins_30000_30500(self):
        fname = "distance_bins_30000_30500"

        p = self._plot_bin_dists(self.d.ld, bin_def="30000 <= distance_bin <= 30500")

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )
        

    # ########## Figure  ##########
    def plot_distance_bins_50000_50500(self):
        fname = "distance_bins_50000_50500"

        p = self._plot_bin_dists(self.d.ld, bin_def="50000 <= distance_bin <= 50500")

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )
        

    # ########## Figure  ##########
    def plot_distance_bins_70000_70500(self):
        fname = "distance_bins_70000_70500"

        p = self._plot_bin_dists(self.d.ld, bin_def="70000 <= distance_bin <= 70500")

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )
        

    # ########## Figure  ##########
    def plot_distance_bins_100000_100500(self):
        fname = "distance_bins_100000_100500"

        p = self._plot_bin_dists(self.d.ld, bin_def="100000 <= distance_bin <= 100500")

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=False
                  )
        

    # ########## Figure  ##########
    def plot_all_contig_len_dist(self):
        fname = "all_contig_len_dist"

        sns.distplot(self.d.contig_info.length, color="coral", kde=0);
        median_contig_len = self.d.contig_info.length.median()
        mean_contig_len = self.d.contig_info.length.mean()
        plt.text(median_contig_len, 400, "*median length = {}".format(median_contig_len), fontsize=14);
        plt.text(mean_contig_len, 200, "*mean length = {}".format(mean_contig_len), fontsize=14);
        plt.text(mean_contig_len, 800, "Number of contigs = {}".format(len(self.d.contig_info)), fontsize=14);
        plt.title("All contigs");

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats
                  )
        plt.close()

    # ### All SNP-pair contigs
    # ########## Figure  ##########
    def plot_all_sp_contig_len_dist(self):
        # ipdb.set_trace()
        fname = "all_sp_contig_len_dist"

        sns.distplot(self.d.sp_contigs.length, color="coral", kde=0);
        sp_median_contig_len = self.d.sp_contigs.length.median()
        sp_mean_contig_len = self.d.sp_contigs.length.mean()
        plt.text(sp_median_contig_len, 100, "*median length = {}".format(sp_median_contig_len), fontsize=14);
        plt.text(sp_mean_contig_len, 50, "*mean length = {}".format(sp_mean_contig_len), fontsize=14);
        plt.text(sp_mean_contig_len, 200, "Number of contigs = {}".format(len(self.d.sp_contigs)), fontsize=14);
        plt.title("Contigs with a SNP-pair");

        self.save_figs(base_dir=self.base_dir,
                       fname=fname,
                       save_types=self.formats
                       )
        plt.close()

    # ########## Figure  ##########
    def plot_distance_VS_contigs_per_bin(self):
        fname = "distance_VS_contigs_per_bin"

        p = self._plot_scat_w_line(gp.ggplot(gp.aes(x='distance_bin', y='contigs_per_bin'), data=self.d.contigs_per_bin))

        # print p

        self.save_figs(base_dir=self.base_dir,
                       fname=fname,
                       save_types=self.formats,
                       is_ggplot=p
                       )
        # 

    # ## All contigs together
    # ########## Figure  ##########
    def plot_distance_VS_r2_all(self):
        # ipdb.set_trace()
        fname = "distance_VS_r2_all"

        p = self._plot_scat_w_line(
            gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=self.d.mean_bin_r2_all.query("distance_bin >= 0")))

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )
        # 

    # ########## Figure  ##########
    def plot_distance_VS_r2_le30K(self):
        fname = "distance_VS_r2_le30K"

        p = self._plot_scat_w_line(
            gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=self.d.mean_bin_r2_all.query("distance_bin <= 30000")))

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )
        # 

    # ########## Figure  ##########
    def plot_distance_VS_r2_le40K(self):
        fname = "distance_VS_r2_le40K"

        p = self._plot_scat_w_line(
            gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=self.d.mean_bin_r2_all.query("distance_bin <= 40000")))

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )
        # 

    # ########## Figure  ##########
    def plot_distance_VS_r2_le50K(self):
        fname = "distance_VS_r2_le50K"

        p = self._plot_scat_w_line(
            gp.ggplot(gp.aes(x='distance_bin', y='R2'), data=self.d.mean_bin_r2_all.query("distance_bin <= 50000")))

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )
        # 

    # ########## Figure  ##########
    def plot_distance_VS_sp_per_bin_gt100(self):
        fname = "distance_VS_sp_per_bin_gt100"

        sp_gt100mask = self.d.len_contigs_per_bin["SNP-pairs"] > 100

        p = self._plot_scat_w_line(gp.ggplot(gp.aes(x='distance_bin',
                                              y='SNP-pairs'),
                                       data=self.d.len_contigs_per_bin[sp_gt100mask])
                             ) + gp.ylab('SNP-pairs per bin')

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )
        # 

    # ########## Figure  ##########
    def plot_distance_VS_sp_per_bin_gt150(self):
        fname = "distance_VS_sp_per_bin_gt150"

        sp_gt150mask = self.d.len_contigs_per_bin["SNP-pairs"] > 150

        p = self._plot_scat_w_line(gp.ggplot(gp.aes(x='distance_bin',
                                              y='SNP-pairs'),
                                       data=self.d.len_contigs_per_bin[sp_gt150mask])
                             ) + gp.ylab('SNP-pairs per bin')

        # print p

        self.save_figs(base_dir=self.base_dir,
                  fname=fname,
                  save_types=self.formats,
                  is_ggplot=p
                  )
        # 

    # In[ ]:
    def plot_distance_VS_avgR2_snpperbin_contigsperbin_q_b_bmin_to_bbmax(self):
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

            p = self._plot_scat_w_line(
                gp.ggplot(
                    gp.aes(
                        x='distance_bin', y='value'),
                    data=self.d.d_bin_v_others_melt.query(
                        "{xmin} <= distance_bin <= {xmax}".format(
                            xmin=xmin,
                            xmax=xmax
                        )
                    )
                )
            ) + \
                gp.facet_wrap("variable") + \
                gp.ggtitle(fname)

            # print p

            self.save_figs(base_dir=self.base_dir,
                           fname=fname,
                           save_types=self.formats,
                           is_ggplot=p
                           )
            # 

    def plot_extras(self):
        # In[ ]:

        # g = sns.PairGrid(d_bin_v_others.loc[:,["R2", "one_minus_cdf_BH",    "SNP-pairs",   "contigs_per_bin",]])
        # g.map_upper(plt.scatter)
        # g.map_lower(sns.kdeplot, cmap="Blues_d")
        # g.map_diag(sns.kdeplot, lw=3, legend=False)
        d_bin_vars = self.d.d_bin_v_others.loc[:, ["R2", "one_minus_cdf_BH", "SNP-pairs", "contigs_per_bin", ]]

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





