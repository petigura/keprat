import numpy as np
import scipy.stats
from matplotlib.pylab import *
import xarray as xr
import keprat.io
import pandas as pd
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

def gaussian(pos, mu, cov):

    """
    Compute the sum of 2D gaussians 

    last dimension of x must equal len(mu)
    """
    assert mu.shape[0]==cov.shape[0]
    assert mu.shape[1]==cov.shape[1]

    out = []
    for i in range(mu.shape[0]):
        _out = scipy.stats.multivariate_normal(mu[i], cov[i]).pdf(pos)
        out.append(_out)
    
    out = np.array(out)
    out = np.sum(out,axis=0)
    return out

class ContourPlotter(object):
    def __init__(self, df):
        self.df = df
        self.dlogP = 0.01
        self.dlogRp = 0.01
        self.Pmin = 1
        self.Pmax = 100
        self.Rpmin = 1
        self.Rpmax = 8
        self.Pticks = [0.3,1,3,10,30,100]
        self.Rpticks = [0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,11.3,16.0]
        self.Pticks = np.array(self.Pticks)
        self.Rpticks = np.array(self.Rpticks)
        
    def compute_density(self):
        logP = arange(log10(self.Pmin),log10(self.Pmax),self.dlogP)
        logRp = arange(log10(self.Rpmin),log10(self.Rpmax),self.dlogRp)
        coords={'logP':logP, 'logRp':logRp}
        _logP,_logRp = meshgrid(logP,logRp,indexing='ij')
        data = {
            'logPc': (['logP', 'logRp'], _logP),
            'logRpc': (['logP', 'logRp'], _logRp),
        }
        ds = xr.Dataset(data,coords=coords)

        n = len(self.df)
        ndim = 2 
        
        # compute midpoints
        mu = np.zeros((n,ndim))
        self.x = log10(self.df.koi_period)
        self.y = log10(self.df.dr25_ror_gdir_srad)
        
        mu[:,0] = self.x
        mu[:,1] = self.y

        # compute uncertainties
        cov = zeros((n,ndim,ndim))
        sigma1 = log10(1.5)
        sigma2 = log10(1 + self.df.dr25_ror_gdir_srad_err1/self.df.dr25_ror_gdir_srad) 
        self.yerr = sigma2
        cov[:,0,0] = sigma1**2
        cov[:,1,1] = sigma2**2

        pos = np.empty(_logP.shape + (2,))
        pos[:, :, 0] = _logP
        pos[:, :, 1] = _logRp
        Z = gaussian(pos, mu, cov)
        ds['Z'] = (['logP','logRp'],Z)
        self.ds = ds

    def plot_contour(self):
        self.ds.Z.plot.contourf(x='logP',cmap=plt.cm.afmhot_r,levels=20)
        xticks(log10(self.Pticks),self.Pticks)
        yticks(log10(self.Rpticks),self.Rpticks)
        xlabel('Orbital Period (days)')
        ylabel('Planet Size (Earth-radii)')
        xlim(log10(self.Pmin),log10(self.Pmax))
        ylim(log10(self.Rpmin),log10(self.Rpmax))
        add_anchored('$N_p$ = {}'.format(len(self.df)),2,prop=dict(size='small'))
        
    def plot_integrated(self):
        pass


    def plot_points(self):
        errorbar(self.x,self.y,yerr=self.yerr,fmt='.',elinewidth=0.5,ms=4)


def plot_before_after(mode):
    import seaborn as sns
    sns.set_context('paper')
    sns.set_style('ticks')
    taucut = ' and (1.1 > tau/tau0 > 0.7)'
    taucutstr = r'$\tau / \tau_0 = 0.7-1.1$'
    df = keprat.io.load_table('cksgaia-planets')
    both = '(dr25_ror_gdir_srad_err1 / dr25_ror_gdir_srad) < 0.1'
    both += ' and gdir_srad <= 10**(0.00025 * (cks_steff - 5500) + 0.2)'
    both += ' and cks_fp==False'
    both += ' and ~(fur17_rcorr_avg > 1.05)'
    both += ' and gaia2_gflux_ratio < 1.1' 
    fig, axL = subplots(ncols=3,nrows=1,figsize=(8,2.5))
    
    if mode=='compare_mult-all_kepmag-lim_mass-all':
        both += ' and kic_kepmag < 14.2'
        both += ' and 0.7 < giso_smass < 1.4'
        _title = '$Kp < 14.2$ and $M_\star = 0.7 - 1.4 M_\odot$'

    elif mode=='compare_mult-all_kepmag-lim_mass-low':
        both += ' and kic_kepmag < 14.2'
        both += ' and 0.7 < giso_smass < 1.0'
        _title = '$Kp < 14.2$ and $M_\star = 0.7 - 1.0 M_\odot$'

    elif mode=='compare_mult-all_kepmag-lim_mass-high':
        both += ' and kic_kepmag < 14.2'
        both += ' and 1.0 < giso_smass < 1.4'
        _title = '$Kp < 14.2$ and $M_\star = 1.0 - 1.4 M_\odot$'

    elif mode=='compare_mult-mult_kepmag-all_mass-all':
        both += ' and 0.7 < giso_smass < 1.4'
        both += ' and multi==True'
        _title = 'Multis and $M_\star = 0.7 - 1.4 M_\odot$'

    elif mode=='compare_mult-mult_kepmag-all_mass-high':
        both += ' and multi==True'
        both += ' and 1.0 < giso_smass < 1.4'
        _title = 'Multis and $M_\star = 1.0 - 1.4 M_\odot$'
    elif mode=='compare_mult-mult_kepmag-all_mass-low':
        both += ' and multi==True'
        both += ' and 0.7 < giso_smass < 1.0'
        _title = 'Multis and $M_\star = 0.7 - 1.0 M_\odot$'

    elif mode=='all':
        mode = [
            'compare_mult-all_kepmag-lim_mass-all',
            'compare_mult-all_kepmag-lim_mass-low',
            'compare_mult-all_kepmag-lim_mass-high',
            'compare_mult-mult_kepmag-all_mass-all',
            'compare_mult-mult_kepmag-all_mass-high',
            'compare_mult-mult_kepmag-all_mass-low',
        ]
        for _mode in mode:
            plot_before_after(_mode)            
    else:
        assert False

    _query1 = both
    _query2 = _query1 + taucut
    samp1 = df.query(_query1)
    samp2 = df.query(_query2)

    sca(axL[0])
    p1 = ContourPlotter(samp1)
    p1.compute_density()
    p1.plot_contour()
    p1.plot_points()

    sca(axL[1])
    p2 = ContourPlotter(samp2)
    p2.compute_density()
    p2.plot_contour()
    add_anchored(taucutstr,3,frameon=False,prop=dict(size='small'))
    title(_title)
    p2.plot_points()

    sca(axL[2])
    temp =  p1.ds['Z'].sum(dim=['logP'])
    plot(temp.logRp,temp / temp.sum())
    xticks(log10(p1.Rpticks),p1.Rpticks)
    temp =  p2.ds['Z'].sum(dim=['logP'])
    plot(temp.logRp,temp /temp.sum())
    xlabel('Planet Size (Earth-radii)')
    tight_layout()
    gcf().savefig('paper/fig_{}.pdf'.format(mode))



def plot_contour(df):
    
    sca(axL[1])
    plot(yi,temp)
    xlabel('log(Rp)')
    ylabel('Density')
    
    sca(axL[2])
    ds['Z'].sum(dim=['logP']).plot()
    xlabel('log(Rp)')
    ylabel('Density (Integrated)')



def band():
    band1 = lambda per : 1.9 * per**-0.05
    band2 = lambda per : 2.3 * per**-0.05
    per = np.logspace(np.log10(0.1),np.log10(300))
    Rp1 = band1(per)
    Rp2 = band2(per)
    fill_between(per,Rp1,Rp2,alpha=0.1,color='c')

def add_anchored(*args,**kwargs):
    ax = gca()
    at = AnchoredText(*args,**kwargs)
    ax.add_artist(at)

def fig_vaneylen_comparison():
    ncols = 4
    nrows = 3
    import seaborn as sns
    sns.set_context('paper',font_scale=1.1)
    sns.set_style('ticks')
    fig,axL = subplots(ncols=ncols,nrows=nrows,figsize=(13,8))
    efmt = dict(fmt='.',elinewidth=0.5)

    df = keprat.io.load_table('cksgaia-planets')
    
    stat1 = [
        [0,0,'v18_Rp','$R_P/R_\star [V18] * R_\star [V18]$'],
        [1,0,'iso_prad','$R_P/R_\star [M15] * R_\star [J17]$'],
        [2,0,'gdir_prad','$R_P/R_\star [M15] * R_\star [F18]$'],
        [0,2,'dr25_ror_gdir_srad','$R_P/R_\star [DR25-MCMC] * R_\star [F18]$'],
        [1,2,'dr25_ror_v18_srad','$R_P/R_\star [DR25-MCMC] * R_\star [V18]$'],
        [2,2,'v18_ror_gdir_srad','$R_P/R_\star [V18] * R_\star [F18]$'],
    ] 
    stat1 = pd.DataFrame(stat1,columns=['irow','icol','y','yann'])
    stat1['xlabel'] = 'Period (days)'
    stat1['ylabel'] = '$R_P$ (Earth-radii)'
    stat1['x'] = 'v18_Period'

    stat2 = [
        [0,1,'1'],
        [1,1,'iso_prad/v18_Rp'],
        [2,1,'gdir_prad/v18_Rp'],
        [0,3,'dr25_ror_gdir_srad/v18_Rp'],
        [1,3,'dr25_ror_v18_srad/v18_Rp'],
        [2,3,'v18_ror_gdir_srad/v18_Rp'],
    ]
    stat2 = pd.DataFrame(stat2,columns=['irow','icol','y'])
    stat2['ylabel'] = '$R_P$ [left] / $R_P$ [V18] '
    stat2['xlabel'] = '$b$ [V18]'
    stat2['x'] = 'v18_b'

    stat = pd.concat([stat1,stat2])
    
    for i, row in stat.iterrows():
        sca(axL[row.irow,row.icol])
        x = df.eval(row.x)
        y = df.eval(row.y)
        setp(gca(),ylabel=row.ylabel,xlabel=row.xlabel)
        errorbar(x,y,**efmt)

    for i, row in stat1.iterrows():
        sca(axL[row.irow,row.icol])
        loglog()
        minorticks_off()
        band()
        xlim(1,100)
        ylim(1,4)
        yt = [1.0,1.4,2.0,2.8,4.0]
        yticks(yt,yt)
        add_anchored(row.yann,loc=2,prop=dict(size='small'))
        

    for i, row in stat2.iterrows():
        sca(axL[row.irow,row.icol])
        y = df.eval(row.y)
        add_anchored(r"RMS = {:.1f}%".format(100*y.std()),loc=2,prop=dict(size='small'))
        ylim(0.8,1.2)

    tight_layout()
    gcf().savefig('paper/fig_v18-comparison-b.pdf')


def plot_compare_ror(cut, yk1, yk2):
    x = 'v18_impact'
    y = yk2 + '/' + yk1
    xerr1 = 'v18_impact_err1'
    xerr2 = 'v18_impact_err2'
    _x = cut.eval(x)
    _y = cut.eval(y)
    _xerr = [-cut.eval(xerr2),cut.eval(xerr1)]
    _label ='N = {}, RMS = {:.1f}%'.format(_y.count(),std(_y)*100)
    errorbar(_x,_y,xerr=_xerr,fmt='o',ms=4,label=_label)

def fig_vanelyen_tau_cuts():
    df = keprat.io.load_table('cksgaia-planets')
    yk1 = 'v18_ror'
    yk2 = 'dr25_RD1_cum'

    def _plot(df):
        plot_compare_ror(df, yk1,yk2)

    fig, axL = subplots(ncols=2,figsize=(8,3.5),sharey=True)
    sca(axL[0])
    title('Multis')
    cut = df.query('multi==True')
    _plot(cut)
    cut = df.query('multi==True and 0.6 < tau/tau0 < 1.1')
    _plot(cut)
    legend()

    sca(axL[1])
    title('Singles')
    cut = df.query('multi==False')
    _plot(cut)
    cut = df.query('multi==False and 0.6 < tau/tau0 < 1.1')
    _plot(cut)
    legend()

    setp(axL[0],ylabel='$(R_p/R_\star)$ [this work] / $(R_p/R_\star)$ [V18]')
    setp(axL,xlabel='$b$ [V18]')
    setp(axL,ylim=(0.8,1.2),xlim=(0,1))
    tight_layout()


def fig_ror_comparison():
    def _plot(cut):
        x = 'v18_impact'
        y = 'dr25_RD1_cum/v18_ror'
        xerr1 = 'v18_impact_err1'
        xerr2 = 'v18_impact_err2'
        _x = cut.eval(x)
        _y = cut.eval(y)
        _xerr = [-cut.eval(xerr2),cut.eval(xerr1)]
        _label ='N = {}, RMS = {:.1f}%'.format(_y.count(),std(_y)*100)
        errorbar(_x,_y,xerr=_xerr,fmt='o',ms=4,label=_label)

    df = keprat.io.load_table('cksgaia-planets')

    fig, axL = subplots()
    sca(axL[0])
