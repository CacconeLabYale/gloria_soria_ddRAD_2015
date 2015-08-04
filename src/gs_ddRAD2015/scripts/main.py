# main.py is part of the 'gloria_soria_ddRAD_2015' package.
# It was written by Gus Dunn and was created on 5/7/15.
#
# Please see the license info in the root folder of this package.

"""
Purpose: provide main entry point to this paper's scripts.

=================================================
main.py
=================================================
"""
import os

import munch
import yaml
import gs_ddRAD2015.scripts.ld_figures
import gs_ddRAD2015.scripts.process_ld_data

__author__ = 'Gus Dunn'

from gs_ddRAD2015 import scripts

import click


def process_config(config):
    """
    Prepare user's config file.

    Also handles validation.

    Args:
        config (str): path to config file.

    Returns:
        conf (dict): configuration values.
    """
    conf = munch.munchify(yaml.load(config))
    conf.meta_data = munch.Munch()
    conf.meta_data.path = os.path.abspath(config.name)

    # # Run all validations that we can do on conf
    # # Can add more later

    return conf

@click.group()
@click.option('--config', default=None,
              help="Path to optional config.yaml",
              type=click.File())
@click.pass_context
def cli(ctx, config):
    """
    Command interface to the paper's LD analysis.

    For command specific help text, call the specific
    command followed by the --help option.
    """

    ctx.obj = munch.Munch()

    if config:
        ctx.obj.CONFIG = process_config(config)



@cli.command()
@click.option('--ld-prog',
              type=click.Choice(['plink', 'vcftools']),
              help="Which program was used to generate the LD values?",
              show_default=True,
              default='vcftools')
@click.option('--distance-bin',
              default=50,
              show_default=True,
              help="How wide do you want the bin window?")
@click.argument('ld_path',
                type=click.Path(exists=True))
@click.argument('out_path',
                type=click.Path())
@click.pass_context
def process_ld_data(ctx, ld_path, out_path, ld_prog, distance_bin):
    """
    Performs the LD analysis.

    \b
    Positional Args:
        ld_path    Path to the table file created by the LD calculation program.
        out_path   Path to where you want to save the results CSV.
    """
    out_dir, out_fname = os.path.split(out_path)

    if not os.path.exists(out_dir):
        click.echo("process_ld_data: Creating directory: {out}.".format(out=out_dir))
        os.makedirs(out_dir)

    gs_ddRAD2015.scripts.process_ld_data.run(ld_path, out_path, ld_prog, distance_bin)




@cli.command()
@click.option('--ld-table',
              type=click.Path(),
              help="path to processed table file.",
              show_default=True,
              default='.')
@click.option('--table-type',
              type=click.Choice(['csv', 'pkl']),
              help="format of table file.",
              show_default=True,
              default='csv')
@click.option('--contig-length',
              type=click.Path(),
              help="path to pickle file.",)
@click.option('--out-dir',
              type=click.Path(dir_okay=True),
              help="path to directory where the figures should go.",
              show_default=True,
              default='.')
@click.option('--formats',
              type=click.Choice(['png', 'svg', 'pdf', 'all', 'none']),
              multiple=True,
              help="the formats that you wish to be produced. Choosing 'none' disables saving of any kind.",
              show_default=True,
              default='none')
@click.option('--save-tables/--no-save-tables',
              show_default=True,
              help="save the data tables in out directory and use them instead of doing the data "
                   "crunching "
                   "over and over each run.",)
@click.option('--force-save/--no-force-save',
              show_default=True,
              help="force overwrite data_tables.", )
@click.pass_context
def ld_figures(ctx, ld_table, table_type, out_dir, formats, contig_length, save_tables, force_save):
    """
    Generates LD figures.

    Takes processed LD d in from of a python pickle and produces figures.
    """

    if not os.path.exists(out_dir):
        click.echo("ld_figures: Creating directory: {out}.".format(out=out_dir))
        os.makedirs(out_dir)

    if 'none' in formats:
        formats = ('none',)

    if 'all' in formats:
        formats = ('png', 'svg', 'pdf')

    gs_ddRAD2015.scripts.ld_figures.run(ld_table, table_type, out_dir, formats, contig_length, save_tables, force_save)


if __name__ == '__main__':
    ld_figures()
