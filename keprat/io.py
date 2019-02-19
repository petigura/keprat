import os
import glob
import cPickle as pickle
from cStringIO import StringIO as sio

import pandas as pd
from astropy.io import ascii
import corner
import numpy as np
from numpy.random import random, multivariate_normal
from chainconsumer import ChainConsumer

DATADIR = os.path.join(os.path.dirname(__file__),'../data/')

def load_table(table, cache=0, cachefn='load_table_cache.hdf', verbose=False):
    """Load tables used in cksmet

    Args:
        table (str): name of table. must be one of
            - nea 


        cache (Optional[int]): whether or not to use the cache
            - 0: don't use the cache recreate all files
            - 1: read from cache
            - 2: write tables to cache

    Returns:
        pandas.DataFrame: table

    """
    if cache==1:
        try:
            df = pd.read_hdf(cachefn,table)
            print "read table {} from {}".format(table,cachefn)
            return df
        except IOError:
            print "Could not find cache file: %s" % cachefn
            print "Building cache..."
            cache=2
        except KeyError:
            print "Cache not built for table: %s" % table
            print "Building cache..."
            cache=2

    if cache==2:
        df = load_table(table, cache=False)
        print "writing table {} to cache".format(table)
        df.to_hdf(cachefn,table,complevel=1,complib='zlib')
        return df

    if table=='kepler-project-mcmc-column-definitions':
        tablefn = os.path.join(DATADIR,'kepler-project-mcmc-column-definitions.txt')
        colspecs = [(0,1),(3,4)]
        df = pd.read_fwf(
            tablefn, comment='#', widths=[20,100],
            names=['column','description']
        )
        
    if table=='coldefs':
        tablefn = os.path.join(DATADIR,'column-definitions.txt')
        colspecs = [(0,1),(3,4)]
        df = pd.read_fwf(
            tablefn, comment='#', widths=[24,100],
            names=['column','description']
        )

    elif table.count('chains'):
        _, dr, id_koicand = table.split('-')
        fmt = {}
        fmt['id_koicand'] = id_koicand
        fmt['id_koi'] = int(id_koicand[1:6])
        fmt['id_plnt'] = int(id_koicand[8:10])
        if dr == 'dr22':
            fn = 'mcmc_chains/dr22/koi{id_koi:}.n/mcmc.{id_koi:}.n{id_plnt:}.dat'.format(**fmt)
            df = read_kepler_tables(fn,'mcmc-dr22')
        elif dr =='dr25':
            fn = 'mcmc_chains/dr25/koi{id_koi:}.n/mcmc.{id_koi:}.n{id_plnt:}.dat.gz'.format(**fmt)
            df = read_kepler_tables(fn,'mcmc-dr25')

    elif table=='v18':
        sing = read_vaneylen('data/vaneylen18/params_table_singles_petigura.txt', mode=2)
        sing = add_prefix(sing,'v18_')
        sing['v18_sample'] = 's'

        mult = read_vaneylen('data/vaneylen15/params_table_multis_petigura.txt', mode=2)
        mult = add_prefix(mult,'v18_')
        mult['v18_sample'] = 'm'
        df = pd.concat([sing,mult])
        
    elif table=='m15':
        import ckscool.io
        df = ckscool.io.load_table('koi-mullally15')
        namemap = {}
        for k in df.columns:
            namemap[k] = k.replace('koi_','m15_')
        df = df.rename(columns=namemap) 

    elif table=='t18':
        import ckscool.io
        df = ckscool.io.load_table('koi-thompson18')
        namemap = {}
        for k in df.columns:
            namemap[k] = k.replace('koi_','t18_')
        df = df.rename(columns=namemap) 

    elif table.count('dr')==1:
        import ckscool.io
        files = glob.glob('mcmc_chains/{}/*.n/mcmc*'.format(table))
        id_koicands = get_id_koicands(files, table)
        mcmc = pd.DataFrame(dict(id_koicand=id_koicands))
        cksc = ckscool.io.load_table('ckscool-targets-cuts')
        cksc = cksc[~cksc.isany]
        cks1 = pd.read_csv('../CKS-Cool/data/cks_physical_merged.csv',index_col=0)
        cols = ['id_koicand']
        cks = pd.concat([cksc[cols],cks1[cols]])
        cks = pd.merge(cks,mcmc)
        cks = cks.sort_values(by='id_koicand').drop_duplicates()
        id_koicands = cks.id_koicand
        #id_koicands = cks.head(30).id_koicand


        df = []
        for id_koicand in id_koicands:
            d = get_summary(id_koicand, table)
            df.append(d)
        df = pd.DataFrame(df)

    elif table=='all':
        dr22 = load_table('dr22', cache=1,cachefn='data/kepler_project_chains.hdf')
        dr25 = load_table('dr25', cache=1,cachefn='data/kepler_project_chains.hdf')
        m15 = load_table('m15', cache=1,)
        t18 = load_table('t18', cache=1)
        v18 = load_table('v18', cache=1)
        df = pd.merge(dr22, dr25, on='id_koicand')
        df = pd.merge(df, m15, on='id_koicand')
        df = pd.merge(df, t18, on='id_koicand')
        df = pd.merge(df, v18, on='id_koicand',how='left')

    elif table=='cksgaia-planets':
        import cksgaia.io
        df = cksgaia.io.load_table(
            'cksgaia-planets',cachefn='../CKS-Gaia/load_table_cache.hdf'
        )
        v18 = load_table('v18',cachefn='load_table_cache.hdf',cache=1)
        dr25 = load_table('dr25',cache=1)
        df = pd.merge(df,v18,on='id_koicand',how='left')
        df = pd.merge(df,dr25,on='id_koicand',how='left')

        df['dr25_ror_v18_srad'] = df.eval('dr25_RD1_cum * v18_srad * 109')
        df['dr25_ror_v18_srad_err1'] = df.eval('dr25_ror_v18_srad * sqrt((dr25_RD1_cum_err1/dr25_RD1_cum)**2 + (v18_srad_err1/v18_srad)**2)')

        df['dr25_ror_gdir_srad'] = df.eval('dr25_RD1_cum * gdir_srad * 109')
        df['dr25_ror_gdir_srad_err1'] = df.eval('dr25_ror_gdir_srad * sqrt((dr25_RD1_cum_err1/dr25_RD1_cum)**2 + (gdir_srad_err1/gdir_srad)**2)')

        df['v18_ror_gdir_srad'] = df.eval('v18_ror * gdir_srad * 109')
        df['v18_ror_gdir_srad_err1'] = df.eval('v18_ror_gdir_srad * sqrt((v18_ror_err1/v18_ror)**2 + (gdir_srad_err1/gdir_srad)**2)')
        df = df.dropna(subset=['gdir_prad'])

        df['tau0'] = 2.036 * df.koi_period**(1/3.) * df.giso_srho**(-1/3.0)
        df['tau'] = df.dr25_TAU1_cum
        df.index = df.id_kic
        df['multi'] = df.groupby('id_kic').size() > 1
        import pdb;pdb.set_trace()
        df = order_columns(df)

    else:
        assert False, "table {} not valid table name".format(table)
    return df

def get_id_koicands(files, dr):
    id_koicands = []
    for f in files:
        try:
            if dr=='dr22':
                _, id_koi, cand, _ = f.split('/')[-1].split('.')
            elif dr=='dr25':
                _, id_koi, cand, _, _ = f.split('/')[-1].split('.')
        except:
            continue

        id_koi = int(id_koi)
        cand = int(cand[1:])
        id_koicand = "K{:05d}.{:02d}".format(id_koi,cand)
        id_koicands.append(id_koicand)

    return id_koicands

from astropy import units as u
from astropy import constants as c

per = 1*u.day
rhostar = 1*u.g/u.cm**3
tau_const = (  
    per**(1/3.) 
    * rhostar**(-1/3.) 
    * 3**(1/3.) 
    * np.pi**(-2/3.) 
    / c.G**(1/3.)
)
tau_const = tau_const.to(u.hour).value 

def get_summary(id_koicand,dr):
    fmt = {}
    fmt['id_koicand'] = id_koicand

    if dr=='dr25':
        nburn = int(4e4)
        table = "chains-dr25-{id_koicand:}".format(**fmt)
    elif dr=='dr22':
        nburn = int(2e4)
        table = "chains-dr22-{id_koicand:}"

    df = load_table(table,cache=1, cachefn='data/kepler_project_chains.hdf')
    df = df.iloc[nburn::10]
    if len(df)==0:
        return {}
     
    df['TAU1'] = tau_const * df.eval(
        '(1 - BB1**2)**(1/2.) * PE1**(1/3.) * RHO**(-1/3.)'
    )

    fgraz = 1.0 * (df['BB1'] > 1.0).sum() / (df['BB1'] > 1.0).count()
    print r"{}% chains with b > 1".format(fgraz * 100)

    d = {}
    
    d['fgraz'] = fgraz
    d['id_koicand'] = id_koicand

    if fgraz==1:
        return d
    df = df[df['BB1'] < 1]
    
    statnames = {'cumulative':'cum'}
    for statistics, stat in statnames.iteritems():
        d3 = {}
        cols = "RHO ZPT EP1 PE1 BB1 RD1 TAU1".split()
        _cols = ["{}_{}_{}".format(dr,k,stat) for k in cols] 
        c = ChainConsumer()
        c.add_chain(np.array(df[cols]), parameters=_cols)
        c.configure(statistics=statistics)
        d2 = c.analysis.get_summary()
        for k in d2.keys():
            lo, mid, hi = d2[k]
            if lo==None:
                lo=mid
            if hi==None:
                hi=mid
            d3[k] = mid
            d3[k+'_err1'] = hi-mid
            d3[k+'_err2'] = lo-mid



        d = dict(d,**d3)
    return d


def read_vaneylen(fn, mode=''):

    if mode==1:
        columns = "0-Kepler, 1-koiname, 2-Ecc, 3-Ecc_low, 4-Ecc_upp, 5-Period, 6-U_Period, 7-Rp, 8-U_Rp, 9-Rovera, 10-Rovera_low, 11-Rovera_upp, 12-RpRs, 13-RpRs_low, 14-RpRs_upp, 15-Mstar, 16-Mstar_low, 17-Mstar_upp, 18-Rstar, 19-Rstar_low, 20-Rstar_upp, 21-Rho, 22-Rho_low, 23-Rho_upp, 24-Temperature, 25-Temperature_low, 26-Temperature_upp, 27-Metl, 28-Metl_low, 29-Metl_upp, 30-Kepmag, 31-relrho, 32-relrho_low, 33-relrho_upp"

    elif mode==2:
        columns = "0-Kepler, 1-koiname, 2-Ecc, 3-Ecc_low, 4-Ecc_upp, 5-Period, 6-U_Period, 7-Rp, 8-U_Rp, 9-Rovera, 10-Rovera_low, 11-Rovera_upp, 12-RpRs, 13-RpRs_low, 14-RpRs_upp, 15-Mstar, 16-Mstar_low, 17-Mstar_upp, 18-Rstar, 19-Rstar_low, 20-Rstar_upp, 21-Rho, 22-Rho_low, 23-Rho_upp, 24-Temperature, 25-Temperature_low, 26-Temperature_upp, 27-Metl, 28-Metl_low, 29-Metl_upp, 30-Kepmag, 31-relrho, 32-relrho_low, 33-relrho_upp, 34-b, 35-b_low, 36-b_upp, 37-ld1, 38-ld1_low, 39-ld1_upp, 40-ld2, 41-ld2_low, 42-ld2_upp"
        
    columns = columns.split(',')
    columns = [c.split('-')[1] for c in columns]
    df = pd.read_csv(fn,sep='\t+',names=columns,skiprows=1,index_col=None)

    for key2 in df.columns:
        if key2.count('low'):
            k = key2.split('_')[0]
            df[k+'_err1'] = df[k+'_upp'] - df[k]
            df[k+'_err2'] = df[k+'_low'] - df[k]
            df = df.drop([k+'_upp',k+'_low'],axis=1)

        if key2.count('U'):
            k = key2.split('_')[1]
            df[k+'_err1'] = df[key2]
            df[k+'_err2'] = - df[key2]
            df = df.drop([key2],axis=1)

    namemap = {
        'Ecc':'ecc','Period':'period','Rp':'prad','Rovera':'rovera','RpRs':'ror','Mstar':'smass','Rstar':'srad','Rho':'srho',
        'Temperature':'steff','Metl':'smet','Kepmag':'kepmag','b':'impact'
    }
    for c in list(df.columns):
        for _old,_new in namemap.iteritems():
            if c==_old or c.count(_old+'_'):
                new = c.replace(_old,_new)
                _d = {c:new}
                df = df.rename(columns=_d)

    df = df.rename(columns={'koiname':'id_koicand'})
    return df

def read_kepler_tables(fn,mode):
    names = load_table('kepler-project-mcmc-column-definitions')
    names = names.column.tolist()
    
    if mode=="mcmc-dr25":
        df = pd.read_table(fn,skiprows=1,sep='\s+',names=names,header=None)
    if mode=="mcmc-dr22":
        df = pd.read_table(fn,skiprows=1,sep='\s+',names=names,header=None)

    return df


def add_prefix(df,prefix,ignore=['id']):
    namemap = {}
    for col in list(df.columns):
        skip=False
        for _ignore in ignore:
            if col.count(_ignore) > 0:
                skip = True
        if not skip:
            namemap[col] = prefix + col 
    df = df.rename(columns=namemap)
    return df

def sub_prefix(df, prefix,ignore=['id']):
    namemap = {}
    for col in list(df.columns):
        skip=False
        for _ignore in ignore:
            if col.count(_ignore) > 0:
                skip = True
        if not skip:
            namemap[col] = col.replace(prefix,'') 
    df = df.rename(columns=namemap)
    return df

def order_columns(df, verbose=False, drop=True):
    columns = list(df.columns)
    coldefs = load_table('coldefs',cache=0)
    cols = []
    for col in coldefs.column:
        if columns.count(col) == 1:
            cols.append(col)

    df = df[cols]
    if verbose and (len(cols) < len(columns)):
        print "table contains columns not defined in coldef"

    return df

