import csv
import itertools
import os
import re

import pandas as pd
import numpy as np
import datetime as dt

import scipy

import munch

import spartan



def recode_villages(df):
    map_func = lambda x: village_id_map[x.upper()]
        
    new_codes = df.Village.apply(map_func)
    df.Village = new_codes

##########################################

def recode_dead(df):
    def recode_func(x):
        # this is treated as an unknown case
        if pd.isnull(x):
            return x

        x = unicode(x)

        # True means DEAD
        # False means LIVE or NOT-DEAD
        # None means unknown

        try:
            # deal with Live type cases
            if x.upper().startswith('L'):
                return False


            if x.startswith('0'):
                return False


            # deal with Dead type cases
            if x.upper().startswith('D'):
                return True


            if x.startswith('1'):
                return True


            # deal with unknown type cases
            if x.upper().startswith('UN'):
                return None
        except AttributeError:
            return x

        msg = "The value {x} was not expected and this function must be corrected to continue.".format(x=x)
        raise ValueError(msg)

    new_dead = df.Dead.apply(recode_func)
    df.Dead = new_dead

##########################################

def recode_teneral(df):
    def recode_func(x):

        # this is treated as an unknown case
        if pd.isnull(x):
            return x

        x = unicode(x)

        # True means teneral
        # False means NOT-teneral
        # None means unknown

        try:
            # deal with NOT-teneral type cases
            if x.upper().startswith('N'):
                return False

            if x.startswith('0'):
                return False

            # deal with Teneral type cases
            if x.upper().startswith('T'):
                return True
                
            if x.startswith('1'):
                return True


            # Deal with unknown type cases
            if x.upper().startswith('UN'):
                return x
        except AttributeError:
            return x

        msg = "The value {x} was not expected and this function must be corrected to continue.".format(x=x)
        raise ValueError(msg)
    
    
    new_teneral = df.Teneral.apply(recode_func)
    df.Teneral = new_teneral

##########################################

def recode_positives(df):
    def recode_func(x):
        # this is treated as an unknown case
        if pd.isnull(x):
            return x

        y = unicode(x)

        # deal with Unknown type cases
        if y.upper().startswith('UN'):
            return None

        if y.upper().startswith('DEAD'):
            return None


        # deal with Positive type cases
        if y.startswith('1'):
            return True


        if y.upper().startswith('TRUE'):
            return True

        if y.upper().startswith('P'):
            return True

        if y.upper().startswith('Y'):
            return True


        # deal with Negative type cases
        if y.upper().startswith('NO'):
            return False

        if y.upper().startswith('FALSE'):
            return False


        if y.startswith('0'):
            return False


        msg = "The value {x} was not expected and this function must be corrected to continue.".format(x=x)
        raise ValueError(msg)


    new_prob = df.prob.apply(recode_func)
    df.prob = new_prob
    
    new_midgut = df.midgut.apply(recode_func)
    df.midgut = new_midgut
    
    new_sal_gland = df.sal_gland.apply(recode_func)
    df.sal_gland = new_sal_gland

##########################################

def recode_species(df):

    recode_func = lambda x: ''.join(x.split('.')).capitalize()

    new_Species = df.Species.apply(recode_func)
    df.Species = new_Species

##########################################

def recode_sex(df):

    recode_func = lambda x: x.upper()

    new_Sex = df.Sex.apply(recode_func)
    df.Sex = new_Sex
    
##########################################

date_delim = re.compile('[\./-]')

def cast_unicode_as_date(x):
    if not isinstance(x, unicode):
        return x
    
    parts = date_delim.split(x)
    
    if len(parts) != 3:
        return x
    
    if len(parts[0]) != 4:
        return x
    
    return dt.datetime(int(parts[0]), int(parts[1]), int(parts[2]))

def recode_date(df):
    new_date = df.Date.apply(cast_unicode_as_date)
    df.Date = new_date

##########################################

fly_no_delim = re.compile('[\W\s]', re.UNICODE)

def split_number(x):
#     ipdb.set_trace()
    
    # to prevent unicode creating a string with a '.' AFTER
    # the numbert we are intersted in!
    try:
        if isinstance(x,float):
            return int(x)
    except ValueError as exc:
        if 'NAN' in exc.message.upper():
            return x
    
    x = unicode(x)
    parts = fly_no_delim.split(x)
    
    try:
        number = int(parts[-1])
        return number
    except ValueError:
        return x


def recode_fly_number(df):
    
    new_fly_number = df.Fly_Number.apply(split_number)
    df.Fly_Number = new_fly_number