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
        out_path   Path to where you want to save the results pickle.
    """
    gs_ddRAD2015.scripts.process_ld_data.run(ld_path, out_path, ld_prog, distance_bin)




@cli.command()
@click.option('--ld-pickle',
              type=click.Path(),
              help="path to pickle file.",
              show_default=True,
              default='.')
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
              help="the formats that you wish to be produced. Choosing 'none' disables saving of any kind.",
              show_default=True,
              default='none')
@click.pass_context
def ld_figures(ctx, ld_pickle, out_dir, formats, contig_length):
    """
    Generates LD figures.

    Takes processed LD d in from of a python pickle and produces figures.
    """

    gs_ddRAD2015.scripts.ld_figures.run(ld_pickle, out_dir, formats, contig_length)


if __name__ == '__main__':
    ld_figures()
