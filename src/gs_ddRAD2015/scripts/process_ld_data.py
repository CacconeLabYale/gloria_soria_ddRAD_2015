import matplotlib.pyplot as plt
import seaborn as sns
import ggplot as gp


import numpy as np
import pandas as pd
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

import click

# import warnings

# warnings.filterwarnings('error')


# this is the command function

echo = click.echo


def run(ld_path, out_path, ld_prog, distance_bin):
    """
    ld_path = "Path to the table file created by the LD calculation program."
    out_path = "Path to where you want to save the results pickle."
    """

    if ld_prog == 'vcftools':
        ld = pd.read_table(ld_path, sep="\t", header=0, names=['CHR_A', 'BP_A', 'BP_B', 'N_INDV', 'R2'])
    else:
        ld = pd.read_table(ld_path, sep=" +", engine='python')

    echo("Calculating distance between snp-pairs.")
    ld['BP_DELTA'] = abs(ld.BP_A - ld.BP_B)

    echo("Establishing the distance_bin column.")
    update_distance_bin(ld, win=distance_bin)

    echo("Re-scaling the R^2 values to (0,1) instead of [0,1].")
    ld["R2_scaled_for_B"] = ld.R2.apply(lambda x: ((x-0.5)*0.999) + 0.5)

    echo("Setting up the models and doing the MAP calculations:")
    models = {}

    for model in yield_models(ld):
        # models[model.bin_id_tag] = model  # Save models in case we need to plot MCMC convergence plots
        record_parameters_and_probabilities(ld, model)

    echo("Saving results to CSV from pandas dataframe.")
    ld.to_csv(path_or_buf=out_path, sep=",", na_rep='nan', float_format=None,
              columns=None, header=True, index=True, index_label=None,
              mode='w', encoding=None, quoting=None,
              quotechar='"', line_terminator='\n', chunksize=None,
              tupleize_cols=False, date_format=None, doublequote=True,
              escapechar=None, decimal='.')

# #################### BELOW ARE HELPER FUNCTIONS #####################


def update_distance_bin(df, win=100):

    assert isinstance(df, pd.DataFrame)

    # add distance_bin column
    # df.loc[:,'distance_bin'] = df.BP_DELTA.apply(lambda x: x - (x % win))
    df['distance_bin'] = df.BP_DELTA.apply(lambda x: x - (x % win))



def yield_models(dataframe):

    bin_ids = dataframe.distance_bin.unique()
    with click.progressbar(bin_ids) as bin_ids:
        for bin_id in bin_ids:

            assert isinstance(bin_id, int)


            # get our distance binned r^2 d in a nice dataframe
            # then drop any rows with R2 == NANs
            data = dataframe.query("(distance_bin == {bin_id})".format(bin_id=bin_id))

            na_mask = data.R2.apply(lambda r2: not np.isnan(r2))
            data = data[na_mask]

            # generate names for stocastics
            alpha_name = "{bin_id}_alpha".format(bin_id=bin_id)
            beta_name = "{bin_id}_beta".format(bin_id=bin_id)
            r2_dist_name = "{bin_id}_r2_distribution_beta".format(bin_id=bin_id)

            # set priors for parameters
            alpha_of_beta = mc.Uniform(alpha_name, 0.01, 10)
            beta_of_beta = mc.Uniform(beta_name, 0.01, 10)

            # set the d
            try:
                r2_distribution_beta = mc.Beta(name=r2_dist_name,
                                               alpha=alpha_of_beta,
                                               beta=beta_of_beta,
                                               value=data.R2_scaled_for_B.dropna(),
                                               observed=True,
                                               verbose=0
                                               )

                # create and yield the model object tagged with its bin_id
                model = mc.Model([r2_distribution_beta, alpha_of_beta, beta_of_beta])
                model.bin_id_tag = bin_id
                model.MCMC_run = None  # allow us to know how many times we had to do the experiments rather than MAP

            except ValueError as exc:
                if "but got (0,)" in exc.message:
                    yield munch.Munch(bin_id_tag=bin_id)

            yield model


def record_parameters_and_probabilities(df, model):

    if not isinstance(model, mc.Model):
        return model

    # set up empty slots for our upcoming values
    # but only if we havent already
    try:
        df.alpha_param[0]
    except AttributeError:
        df['alpha_param'] = np.nan
        df['beta_param'] = np.nan
        df['cdf'] = np.nan
        df['one_minus_cdf'] = np.nan
        df['one_minus_cdf_BH'] = np.nan
        df['MAP_succeeded'] = 'no val'

    # set up a mask for this distance_bin
    bin_mask = df.distance_bin == model.bin_id_tag

    # perform Max A Posteriori
    model_runner = mc.MAP(model)


    # Parameter names for access through the model
    a_name = "{bin_id}_alpha".format(bin_id=model.bin_id_tag)
    b_name = "{bin_id}_beta".format(bin_id=model.bin_id_tag)

    try:
        # if we can take the short cut, great!
        model_runner.fit()
        model.MCMC_run = False
        df.loc[bin_mask, 'MAP_succeeded'] = "MAP"
    except (RuntimeError, ZeroDivisionError) as exc:

        run_mcmc = False
        dbin = df.query(""" distance_bin == {d} """.format(d=model.bin_id_tag))

        if "Posterior probability optimization converged to value with zero probability." in exc.message:
            # Otherwise lets run the MCMC experiment
            echo("\nEncountered 'zero probability' error, running MCMC for distance_bin: {bin_id} with {num} SNPs.".format(bin_id=model.bin_id_tag,num=len(dbin)))

            run_mcmc = True

        elif "float division by zero" in exc.message:
            # Otherwise lets run the MCMC experiment
            echo("\nEncountered 'float division by zero' error, running MCMC for distance_bin: {bin_id} with {num} SNPs.".format(bin_id=model.bin_id_tag,num=len(dbin)))

            run_mcmc = True

        else:
            raise

        if run_mcmc:
            run_MCMC_if_MAP_fails(model,a_name,b_name)
            df.loc[bin_mask, 'MAP_succeeded'] = "MCMC"


    # get and store parameters
    alpha, beta = model.get_node(a_name),model.get_node(b_name)

    df.loc[bin_mask,"alpha_param"] = alpha.value.item()
    df.loc[bin_mask,"beta_param"] = beta.value.item()

    # get and store probabilities
    # Record the probabilities of obtaining each R2 value (or smaller) in the bin afer scaling
    df.loc[bin_mask, 'cdf'] = scipy.stats.beta.cdf(df[bin_mask].R2_scaled_for_B,
                                                   alpha.value.item(),
                                                   beta.value.item())
    df.loc[bin_mask, 'one_minus_cdf'] = 1 - df[bin_mask]['cdf']
    df.loc[bin_mask, 'one_minus_cdf_BH'] = smm.multipletests(df[bin_mask].one_minus_cdf,
                                                         method='fdr_bh')[1]


def run_MCMC_if_MAP_fails(model,a_name,b_name):
    model_runner = mc.MCMC(model)
    model.MCMC_run = model_runner

    # re-initialize the alpha and beta starting values randomly
    # to avoid non-allowed values that seem to slip in from the MAP.fit()
    current_alpha,current_beta = model.get_node(a_name),model.get_node(b_name)
    current_alpha.random()
    current_beta.random()

    # do the learning
    model_runner.sample(iter=40000, burn=20000, thin=1)

if __name__ == '__main__':
    run()
