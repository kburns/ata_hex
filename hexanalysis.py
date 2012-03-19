"""
Plotting and analysis procedures.

Author: Keaton J. Burns <keaton.burns@gmail.com>
Affiliation: UC Berkeley

"""


import os
import numpy as np
import scipy
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
import hextoolkit
from kbtoolkit import normal_dist


def abs_dist(tag, percentile=95, flag=False, plot=True, saveas=False):
    """
    Determine cutoff percentile for specified tag, plot histogram, flag
    
    Inputs:
        tag             tag for data to examine
        percentile      percentile for cutoff, in range [0, 100]
        flag            flagging number if flagging obs past cutoff is desired
        plot            Boolean to create histogram plot
        save            Boolean to save histogram plot

    """
    
    
    # Get tag data
    data = hextoolkit.querysql([tag], exclude_flagged=False)
    tagdata = abs(data[tag])
    
    # Use scipy percentile function to get proper percentile
    cutoff = scipy.percentile(tagdata, percentile)
    bins = 100
    
    if plot:
        # Make histogram of Chi-Squared Values
        fig = plt.figure(1, figsize=(8,7))
        fig.clf()
        
        # Full histogram
        ax1 = fig.add_subplot(211)
        ax1.hist(tagdata, bins=bins, log=True)
        ax1.set_title('Distribution of ' + tag + ' magnitude')
        ax1.set_xlabel(tag)
        ax1.set_ylabel('Count')
        ax1.set_ylim([0.5,ax1.axis()[3]])
        
        # Zoom near cutoff
        ax2 = fig.add_subplot(212)
        ax2.hist(tagdata, bins=bins, log=True, range=(0, cutoff * 2))
        ax2.set_xlabel(tag)
        ax2.set_ylabel('Count')
        ax2.axvline(cutoff, color='r', label=' %i percentile: %f' %(percentile, cutoff))
        ax2.set_ylim([0.5,ax2.axis()[3]])
        
        # Small font legend indicated cutoff
        prop = matplotlib.font_manager.FontProperties(size='small')
        plt.legend(prop=prop)
        
        if saveas:
            if saveas is True: saveas = 'dist_' + tag + '.png'
            plt.savefig(saveas)
    
    if flag:
        # Default flag value
        if flag is True:
            flag = 10
            
        hextoolkit.flagsql('WHERE oid IN (SELECT oid FROM runs natural join obs WHERE abs(' + tag + ') >= %f)' %cutoff, num=flag)
     
 
def flagging(analyze=False, saveas=False):
    """Store flag settings"""
    
    # Reset all flags
    hextoolkit.flagsql(num=0, add=False)
    
    # Dont use 700
    hextoolkit.flagsql(wherecmd="WHERE rid IN (SELECT rid FROM runs WHERE freq = 700)", num=2**15)
    
    flaglist = [('squintaz*freq', 99, 2**0),
                ('squintaz_uc*freq', 98, 2**1),
                ('squintel*freq', 99, 2**2),
                ('squintel_uc*freq', 98, 2**3),
                ('sefd', 98, 2**4),
                ('sumchisq', 90, 2**5),
                ('x_width_el*freq', 99, 2**6),
                ('x_width_el_uc*freq', 98, 2**7),
                ('x_width_az*freq', 99, 2**8),
                ('x_width_az_uc*freq', 98, 2**9),
                ('y_width_el*freq', 99, 2**10),
                ('y_width_el_uc*freq', 98, 2**11),
                ('y_width_az*freq', 99, 2**12),
                ('y_width_az_uc*freq', 98, 2**13)]
                
    for t,p,f in flaglist:
        abs_dist(t, p, flag=f, plot=False)
        
    if analyze:
        if analyze is True: analyze = 'freq'
        data = hextoolkit.querysql([analyze, 'flag'], exclude_flagged=False)
        
        a_list = np.sort(np.unique(data[analyze]))
        N = np.size(a_list)
        M = len(flaglist)
        flaggrid = np.zeros((N,M+1))
        flag_totals = np.zeros(N)
        a_count = np.zeros(N)
        
        for row in data:
             for n in xrange(N):
                if row[analyze] == a_list[n]:
                    a_count[n] += 1
                    if row['flag'] != 0:
                        flag_totals[n] += 1
                        for m in xrange(M):
                            if np.bitwise_and(row['flag'], flaglist[m][2]) != 0:
                                flaggrid[n, m] += 1
                    continue
                            
        flaggrid[:, -1] = flag_totals
        
        print a_list
        print 'Totals:', flag_totals / np.asfarray(a_count)
                    
        plt.figure(figsize=(12,8))
        plt.clf()
        ax = plt.gca()
        
        flaggrid = (flaggrid.T/np.asfarray(a_count)).T
        plt.imshow(flaggrid, interpolation='nearest')
        
        for n in xrange(N):
            for m in xrange(M + 1):
                plt.text(m, n, '%.2f' %flaggrid[n, m], fontsize=8, color='w', horizontalalignment='center', verticalalignment='center', weight='bold')
        
        plt.setp(ax, xticks=np.arange(M+1))
        xticklabs = [row[0] for row in flaglist]
        xticklabs.append('total')
        xtickNames = plt.setp(ax, xticklabels=xticklabs)
        plt.setp(xtickNames, rotation=80, fontsize=8)
        plt.xlabel('Flag')
        
        plt.setp(ax, yticks=np.arange(N))
        ytickNames = plt.setp(ax, yticklabels=a_list)
        plt.setp(ytickNames, fontsize=8)
        plt.ylabel(analyze)
        
        cbar = plt.colorbar()
        cbar.set_label('Fraction Flagged')
        plt.draw()
        
        if saveas:
            if saveas is True: saveas = 'current_flagging.png'
            plt.savefig(saveas)
        
    
def stats():
    """Print basic statistics"""

    data = hextoolkit.querysql(['antnum', 'feed', 'feedrev', 'freq', 'flag'], 
                               exclude_flagged=False)
                               
    antlist = np.unique(data['antnum'])
    feedlist = np.unique(data['feed'])
    
    tot_obs = data.size
    tot_flagged = np.sum(data['flag'] != 0)
    print 'Total Observations = ', tot_obs
    print 'Flagged Observations = ', tot_flagged, np.float(tot_flagged) / tot_obs
    print 'Num Antennas = ', antlist.size
    print 'Num Feeds = ', feedlist.size
    
    # Find antennas with multiple feeds
    ant_mult_feeds = 0
    amf_obs = 0
    for ant in antlist:
        adata = data[data['antnum'] == ant]
        if np.unique(adata['feedrev']).size >= 2:
            amf_obs += adata.size
            ant_mult_feeds += 1
    
    print 'Num Antennas with Multiple Feeds = ', ant_mult_feeds
    print '    Obs                          = ', amf_obs
    
    # Find feeds on multiple antennas
    feeds_mult_ant = 0
    fma_obs = 0
    for feed in feedlist:
        fdata = data[data['feed'] == feed]
        if np.unique(fdata['antnum']).size >= 2:
            fma_obs += fdata.size
            feeds_mult_ant += 1
    
    print 'Num Feeds with Multiple Antennas = ', feeds_mult_ant
    print '    Obs                          = ', fma_obs


def beam_shape_plot(coords='mag_angle', type='scatter', saveas=False, wherecmd=''):
    """
    Make distribution plots of beam shape.
    
    Inputs:
        coords      One of ['az_el', 'mag_angle', 'mag_ecc']
        type        One of ['scatter','density']
        saveas      Path and filename for saving (will be prepended with polarization)
        wherecmd    "Where" filter for data
    
    """
    
    if coords == 'az_el':
        data = hextoolkit.querysql(['freq', 'x_width_el', 'x_width_az',
                                    'y_width_el', 'y_width_az'], wherecmd=wherecmd)
        fghz = data['freq'] / 1000.
        
        xx = data['x_width_az'] * 2 * fghz
        xx_uc = data['x_width_az_uc'] * 2 * fghz
        xy = data['x_width_el'] * 2 * fghz
        xy_uc = data['x_width_el_uc'] * 2 * fghz
        yx = data['y_width_az'] * 2 * fghz
        yx_uc = data['y_width_az_uc'] * 2 * fghz
        yy = data['y_width_el'] * 2 * fghz
        yy_uc = data['y_width_el_uc'] * 2 * fghz
        
        hline = 3.5
        vline = 3.5
        
        xlabel = lambda pol: r'$\rm{%s \, Azimuth \, FWHM \, (degrees \times \, f_{GHz})}$' %pol
        ylabel = lambda pol: r'$\rm{%s \, Elevation \, FWHM \, (degrees \times \, f_{GHz})}$' %pol
        
        if saveas is True:
            saveas = 'az_el.png'
            
        xlim = None
        ylim = None
        
        print 'DATA: xwa, ywa, xwe, ywe'
    
    elif coords == 'mag_angle':
        data = hextoolkit.querysql(['freq', 'xmag', 'xangle', 'ymag', 'yangle'], wherecmd=wherecmd)
        fghz = data['freq'] / 1000.
        
        xx = data['xmag'] * 2 * fghz
        xx_uc = data['xmag_uc'] * 2 * fghz
        xy = data['xangle'] * 180. / np.pi
        xy_uc = data['xangle_uc'] * 180. / np.pi
        yx = data['ymag'] * 2 * fghz
        yx_uc = data['ymag_uc'] * 2 * fghz
        yy = data['yangle'] * 180. / np.pi
        yy_uc = data['yangle_uc'] * 180. / np.pi
        
        hline = 45.0
        vline = 3.5
        
        xlabel = lambda pol: r'$\rm{%s \, Beam \, Size \, (degrees \times \, f_{GHz})}$' %pol
        ylabel = lambda pol: r'$\rm{%s \, Beam \, Angle \, (degrees)}$' %pol
        
        if saveas is True:
            saveas = 'mag_angle.png'
            
        xlim = None
        ylim = None
        
        print 'DATA: x mag, y mag, x angle, y angle'
    
    elif coords == 'mag_ecc':
        data = hextoolkit.querysql(['freq', 'xmag', 'ymag'], wherecmd=wherecmd)
        fghz = data['freq'] / 1000.
        
        xm = data['xmag']
        xm_uc = data['xmag_uc']
        
        xa2 = data['x_width_az'] ** 2
        xa2_uc = 2 * data['x_width_az'] * data['x_width_az_uc']
        xe2 = data['x_width_el'] ** 2
        xe2_uc = 2 * data['x_width_el'] * data['x_width_el_uc']
        
        pmask = (xa2 >= xe2)
        xfrac = xa2 / xe2
        xfrac[pmask] = (xe2 / xa2)[pmask]
        xfrac_uc = xfrac * np.sqrt((xa2_uc / xa2)**2 + (xe2_uc / xe2)**2)
        xecc = np.sqrt(1 - xfrac)
        xecc[pmask] *= -1
        xecc_uc = 0.5 * xecc / (1 - xfrac) * xfrac_uc
        
        ym = data['ymag']
        ym_uc = data['ymag_uc']
        
        ya2 = data['y_width_az'] ** 2
        ya2_uc = 2 * data['y_width_az'] * data['y_width_az_uc']
        ye2 = data['y_width_el'] ** 2
        ye2_uc = 2 * data['y_width_el'] * data['y_width_el_uc']
        
        pmask = (ya2 >= ye2)
        yfrac = ya2 / ye2
        yfrac[pmask] = (ye2 / ya2)[pmask]
        yfrac_uc = yfrac * np.sqrt((ya2_uc / ya2)**2 + (ye2_uc / ye2)**2)
        yecc = np.sqrt(1 - yfrac)
        yecc[pmask] *= -1
        yecc_uc = 0.5 * yecc / (1 - yfrac) * yfrac_uc
        
        xx = data['xmag'] * 2 * fghz
        xx_uc = data['xmag_uc'] * 2 * fghz
        xy = xecc
        xy_uc = xecc_uc
        yx = data['ymag'] * 2 * fghz
        yx_uc = data['ymag_uc'] * 2 * fghz
        yy = yecc
        yy_uc = yecc_uc
        
        hline = 0.0
        vline = 3.5
        
        xlabel = lambda pol: r'$\rm{%s \, Beam \, Size \, (degrees \times \, f_{GHz})}$' %pol
        ylabel = lambda pol: r'$\rm{%s \, Eccentricity}$' %pol
        
        if saveas is True:
            saveas = 'mag_ecc.png'
            
        xlim = None
        ylim = [-1, 1]
        
        print 'DATA: x mag, y mag, x ecc, y ecc'

    # Print stats
    print 'MEDIAN:', np.median(xx), np.median(yx), np.median(xy), np.median(yy)
    print 'MEAN:  ', np.mean(xx), np.mean(yx), np.mean(xy), np.mean(yy)
    print ' ERR:  ', np.std(xx) / np.sqrt(xx.size), np.std(yx) / np.sqrt(yx.size), np.std(xy) / np.sqrt(xy.size), np.std(yy) / np.sqrt(yy.size)
    print 'STDEV: ', np.std(xx), np.std(yx), np.std(xy), np.std(yy)
    
    # Definitions for the axes
    nullfmt = matplotlib.ticker.NullFormatter()
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    for pol in ['X', 'Y']:

        plt.figure(1, figsize=(8,8))
        plt.clf()
        
        if pol == 'X':
            x, xuc, y, yuc = xx, xx_uc, xy, xy_uc
        else:
            x, xuc, y, yuc = yx, yx_uc, yy, yy_uc
                
        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        
        # Remove histogram ticklabels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        # Scatter/Density plot
        if type == 'scatter':
            axScatter.errorbar(x, y, xerr=xuc, yerr=yuc, fmt='o')
        if type == 'density':
            density, xedges, yedges = np.histogram2d(x, y, bins=100)
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            logd = np.log10(density)
            logd[logd == -np.inf] = -1
            levels = np.arange(-0.5, np.ceil(np.max(logd)), 0.5)
            axScatter.contour(logd.T, levels, extent=extent)
         
        # Histograms
        bins = 40
        xn, xb, xp = axHistx.hist(x, bins=bins, log=True)
        yn, yb, yp = axHisty.hist(y, bins=bins, orientation='horizontal', log=True)
                
        # Limits
        if xlim:
            axScatter.set_xlim(xlim)
        else:
            xlim = axScatter.get_xlim()
        if ylim:
            axScatter.set_ylim(ylim)
        else:
            ylim = axScatter.get_ylim()
          
        axHistx.set_xlim(xlim)
        axHisty.set_ylim(ylim)
        
        # Gaussian fits
        N = 1000
        gx = np.linspace(xlim[0], xlim[1], N)
        gy = np.linspace(ylim[0], ylim[1], N)
        xylim = axHistx.get_ylim()
        yxlim = axHisty.get_xlim()
        fx = normal_dist(np.mean(x), np.std(x))(gx) * x.size * np.diff(xb)[0]
        fy = normal_dist(np.mean(y), np.std(y))(gy) * y.size * np.diff(yb)[0]
        axHistx.plot(gx, fx, 'g')
        axHisty.plot(fy, gy, 'g')
        axHistx.set_ylim(xylim)
        axHisty.set_xlim(yxlim)
        
        # Axis labels
        axScatter.set_xlabel(xlabel(pol))
        axScatter.set_ylabel(ylabel(pol))
        
        # Lines
        axScatter.axhline(hline, c='r')
        axScatter.axvline(vline, c='r')
        axHistx.axvline(vline, c='r')
        axHisty.axhline(hline, c='r')

        # Save
        if saveas:
            folder, file = os.path.split(saveas)
            plt.savefig(os.path.join(folder, pol + '_' + file))


def beam_shape_by_freq(coords='az_el', type='density', saveas=False):
    """Make beam shape plots for each frequency."""

    data = hextoolkit.querysql(['freq'])
    freqs = np.unique(data['freq'])
    
    if saveas is True:
        saveas = 'beam_shape.png'
    
    for f in freqs:
        fsave = str(f) + '_' + saveas
        beam_shape_plot(coords=coords, type=type, saveas=fsave, wherecmd='WHERE freq=%f' %f)


def squint_stddev(rev=False):
    
    # Query data
    data = hextoolkit.querysql(['antfeedrev','date','squintmag','squintangle','freq'])
    
    antfeedrevs = np.unique(data['antfeedrev'])
    freqlist = np.unique(data['freq'])
    daylist = np.unique(np.round(data['date']))

    magdev = []
    angledev = []
    
    # Examine each antfeedred
    for afr in antfeedrevs:
        if rev:
            if np.round(np.mod(afr, 1), 1) != rev: continue

        idata = data[np.where(data['antfeedrev'] == afr)]

        # Examine each day
        for f in freqlist:
            ifdata = idata[np.where(idata['freq'] == f)]
            # Skip if no points for this antfeedrev in this run
            if ifdata.size < 2: 
                continue
            
            # Compute magnitude standard deviation as a fraction of the mean
            mag = ifdata['squintmag']
            mag /= mag.mean()
            magdev.append(np.std(mag))
            
            # Compute angle standard deviation as a fraction of 2*pi, minimizing over circle
            angle = ifdata['squintangle']
            modangles = np.linspace(0, 2 * np.pi, 100)
            moddevs = []
            for i in modangles:
                moddevs.append(np.std(np.mod(angle + i, 2 * np.pi)))
            angledev.append(np.min(moddevs) / (np.pi / 2))

    return np.array(magdev), np.array(angledev)
    
    
def squint_stddev_plots(bins=40, save=False):

    rev = [0.1, 0.2, 0.3]

    # Plot histograms of fits
    fig = plt.figure(1, figsize=(10,10))
    fig.clf()
    bins = np.linspace(0, 2, bins)

    for i,r in enumerate(rev):
        magdev, angledev = squint_time_stddev(rev=r)
                        
        # Power law histogram
        ax = fig.add_subplot(3, 1, i + 1)
        #ax.hist([magdev, angledev], bins=bins, histtype='step', color=['b', 'r'],
        #        label=[r'$\rm{Mag./Mean}$', r'$\rm{Angle/(\pi / 2)}$'], linewidth=2.0)
        ax.hist(magdev, bins=bins, histtype='step', color='b', label=r'$\rm{Mag./Mean}$', lw=2.0)
        ax.hist(angledev, bins=bins, histtype='step', color='r', label=r'$\rm{Angle/(\pi / 2)}$', lw=2.0, ls='dashed')
        ax.set_ylabel(r'$\rm{Count}$')
        ax.text(0.75, 0.8, r'$\rm{Feed Rev} = %.1f$' %r, transform=ax.transAxes)
        ax.text(0.75, 0.7, r'$\rm{Mag. \ Mean} = %.2f \pm %.2f$' %(np.mean(magdev), np.std(magdev) / np.sqrt(magdev.size)), transform=ax.transAxes)
        ax.text(0.75, 0.6, r'$\rm{Angle \ Mean} = %.2f \pm %.2f$' %(np.mean(angledev), np.std(angledev) / np.sqrt(angledev.size)), transform=ax.transAxes)

        if i == 2:
            ax.set_xlabel(r'$\rm{Normalized \ Standard \ Devation}$')
            ax.legend(loc='lower right')
            
    if save:
        if save is True: save = 'squinttime_rev.png'
        plt.savefig(save)




## May require clean up

def magfreq_plot(mincount=1, log=False, saveas='magfreq.pdf'):
    """Study magnitude vs frequency"""
    
    # Setup pdf output
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(saveas)
    
    # Query data
    data = hextoolkit.querysql(['antfeedrev','date','squintmag','freq'])
    
    antfeedrevs = np.unique(data['antfeedrev'])
    freqlist = np.unique(data['freq'])
    
    # Examine each antfeedred
    for afr in antfeedrevs:
        plt.figure(1, figsize=(9, 9))
        plt.clf()
        
        idata = data[np.where(data['antfeedrev'] == afr)]
        
        # Examine each frequency
        for f in freqlist:
            ifdata = idata[np.where(idata['freq'] == f)]
            # Skip if not enough points for this frequency
            if ifdata.size < mincount: continue
            
            # Get statistics of data closest to median
            ifmedianmag = np.median(ifdata['squintmag'])
            ifindex = (np.abs(ifdata['squintmag'] - ifmedianmag)).argmin()
            ifmag = ifdata['squintmag'][ifindex]
            iferr = ifdata['squintmag_uc'][ifindex]
            
            # Modify errors to plot on log scale (low bar must be positive)
            if log: 
                plt.yscale('log')
                plt.xscale('log')
                iferr = np.array([[min(iferr, ifmag - 1e-5)], [iferr]])

            plt.errorbar(f, ifmag, yerr=iferr, fmt='o', label=str(ifdata.size))
    
        plt.title('antfeedrev = ' + str(afr))
        plt.xlabel('freq')
        plt.ylabel('squintmag')
        
        plt.axis([500, 8e3, 0.01, 35])
    
        # Small legend
        prop = matplotlib.font_manager.FontProperties(size='small')
        plt.legend(prop=prop, numpoints=1, title='Observations')
    
        pp.savefig()
        
    pp.close()
  
  
def magfreq_powerlaw_fit(mincount=2, rev=False):
    """Compute power laws of squint vs mag dependence for each antfeedrev"""

    # Query data
    data = hextoolkit.querysql(['antfeedrev','date','squintmag','freq'])
    
    antfeedrevs = np.unique(data['antfeedrev'])
    freqlist = np.unique(data['freq'])
    daylist = np.unique(np.round(data['date']))

    powers = []
    chisqbest = []
    chisq0 = []
    chisq1 = []
    sizes = []

    # Examine each antfeedred
    for afr in antfeedrevs:
        if rev:
            if np.round(np.mod(afr, 1), 1) != rev: continue

        idata = data[np.where(data['antfeedrev'] == afr)]
        
        # Examine each day
        for d in daylist:
            iddata = idata[np.where(np.round(idata['date']) == d)]
            # Skip if no points for this antfeedrev in this run
            if iddata.size < mincount: 
                continue
            else:
                sizes.append(iddata.size)
            
            idf = iddata['freq']
            ids = iddata['squintmag']
            ids_uc = iddata['squintmag_uc']
            
            #if iddata.size > 10:
            #    plt.clf()
            #    plt.errorbar(idf, ids, yerr=ids_uc, fmt='.')
            #    raw_input('Paused')
            
            # Compute best fit power law index   
            logfit = np.polyfit(np.log10(iddata['freq']), 
                                np.log10(iddata['squintmag']), 1)
            powers.append(logfit[0])                    
            fit = 10**logfit[1] * idf ** logfit[0]
            chisq = np.sum((ids - fit)**2 / ids_uc**2)
            chisqbest.append(chisq)
            
            # Test 0 power law fit
            fit0 = np.mean(ids)
            chisq = np.sum((ids - fit0)**2 / ids_uc**2)
            chisq0.append(chisq)
            
            # Test -1 power law fit
            A = np.mean(ids * idf)
            fit1 = A * 1./ idf
            chisq = np.sum((ids - fit1)**2 / ids_uc**2)
            chisq1.append(chisq)
 
    return np.array(powers), np.array(chisqbest), np.array(chisq0), np.array(chisq1), np.array(sizes)
    
        
def magfreq_powerlaw_rev(mincount=2, bins=50, save=False):
    """Compute power laws by rev of squint vs mag dependence for each antfeedrev
    
    -look at large amplitudes
    -test hypothesis -1 fit
    -goodness of fit
    
    """

    rev = [0.1, 0.2, 0.3]

    # Plot histograms of fits
    fig = plt.figure(1, figsize=(10,10))
    fig.clf()

    for i in xrange(len(rev)):
        powers, chb, ch0, ch1, sizes = magfreq_powerlaw_fit(mincount=mincount, rev=rev[i])
        
        # Disregard outliers
        #outliers = (chb / sizes >= 10)
        #outliers = (np.abs(powers) >= 5)
        outliers = (ch0 >= 100) + (ch1 >= 100)
        print '%i outliers of %i' %(outliers.sum(), powers.size)
        powers = powers[~outliers]
        ch0 = ch0[~outliers]
        ch1 = ch1[~outliers]
        sizes = sizes[~outliers]
        
        # Print likelihood ratios
        pdf0 = scipy.stats.distributions.chi2.pdf(ch0, sizes)
        pdf1 = scipy.stats.distributions.chi2.pdf(ch1, sizes)
        print 'Log-Likelihood ratio (1/0):', np.log10(np.product(pdf1/pdf0))
                        
        # Power law histogram
        ax = fig.add_subplot(3, 1, i + 1)
        ax.hist(powers, bins=np.linspace(-5, 5, bins))
        title_str = r'count $\geq$ %i, rev = %.1f' %(mincount, rev[i])
        #ax.set_title(title_str)
        ax.set_ylabel(r'$\rm{Count}$')
        #ax.set_ylim([0.5,ax.axis()[3]])
        #ax.axis([-2.5,0.5,0,6])
        ax.text(0.75, 0.8, r'$\rm{Feed Rev} = %.1f$' %rev[i], transform=ax.transAxes)
        ax.text(0.05, 0.8, r'$\rm{Mean} = %.2f \pm %.2f$' %(np.mean(powers), np.std(powers) / np.sqrt(powers.size)), transform=ax.transAxes)
        ax.text(0.05, 0.7, r'$\rm{StDev} = %.2f$' %np.std(powers), transform=ax.transAxes)

        if i == 2:
            ax.set_xlabel(r'$\rm{Power \ law \ index} \ \alpha \ (|\vec{S}| \propto f^{\alpha})$')
            
    if save:
        if save == True: save = 'powerlaw_rev.png'
        plt.savefig(save, type='png')
 
 
def ant_corr():
    """Compute correlations between median squintmag for each frequency and ant/feed number"""
    
    # Query data
    data = hextoolkit.querysql(['antnum', 'feed', 'squintmag', 'freq'])
    
    antlist = np.unique(data['antnum'])
    feedlist = np.unique(data['feed'])
    freqlist = np.unique(data['freq'])
    
    ant_corr = []
    
    fig = plt.figure(1, figsize=(16,12))
    plt.clf()
    A = 0
    
    for ant in antlist:
        
        adata = data[np.where(data['antnum'] == ant)]
        adict={'freq':[], 'amp':[]}
        if np.unique(adata['feed']).size < 2: continue

        A += 1
        ax = fig.add_subplot(3,3,A)

        for f in np.unique(adata['feed']):
            afdata = adata[np.where(adata['feed'] == f)]
            adict[str(f) + '_freq'] = []
            adict[str(f) + '_amp'] = []
            
            for fr in np.unique(afdata['freq']):
                affdata = afdata[np.where(afdata['freq'] == fr)]
                adict[str(f) + '_freq'].append(fr)
                adict[str(f) + '_amp'].append(np.median(affdata['squintmag']))
                adict['freq'].append(fr)
                adict['amp'].append(np.median(affdata['squintmag']))
               
            ax.loglog(adict[str(f) + '_freq'], adict[str(f) + '_amp'], '.-')
            if A == 7:
                ax.set_xlabel('Freq')
                ax.set_ylabel('Sq Mag')
        
        pearsonr = scipy.stats.pearsonr(np.log10(adict['freq']), np.log10(adict['amp']))
        ant_corr.append(pearsonr[0])
        
        ax.set_title('Ant %i, Corr = %f' %(ant, pearsonr[0]))
        
    plt.draw()
    plt.savefig('ant_feed_corr/ant_mult_feeds.png', type='png')
        
    fig = plt.figure(2)
    plt.clf()
    plt.hist(ant_corr, bins=20)
    plt.text(-0.3, 1.5, 'Mean = %f' %np.average(ant_corr))
    plt.axis([-1, 0, 0, 2])
    plt.savefig('ant_feed_corr/ant_corr.png', type='png')
 
 
def feed_corr():
    """Compute correlations between median squintmag for each frequency and ant/feed number"""
    
    # Query data
    data = hextoolkit.querysql(['antnum', 'feed', 'squintmag', 'freq'])
    
    antlist = np.unique(data['antnum'])
    feedlist = np.unique(data['feed'])
    freqlist = np.unique(data['freq'])
    
    ant_corr = []
    
    fig = plt.figure(1, figsize=(16,12))
    plt.clf()
    A = 0
    
    for feed in feedlist:
        adata = data[np.where(data['feed'] == feed)]
        adict={'freq':[], 'amp':[]}
        if np.unique(adata['antnum']).size < 2: continue
        
        A += 1
        ax = fig.add_subplot(3,3,A)

        for ant in np.unique(adata['antnum']):
            afdata = adata[np.where(adata['antnum'] == ant)]
            adict[str(ant) + '_freq'] = []
            adict[str(ant) + '_amp'] = []
            
            for fr in np.unique(afdata['freq']):
                affdata = afdata[np.where(afdata['freq'] == fr)]
                adict[str(ant) + '_freq'].append(fr)
                adict[str(ant) + '_amp'].append(np.median(affdata['squintmag']))
                adict['freq'].append(fr)
                adict['amp'].append(np.median(affdata['squintmag']))
                
            ax.loglog(adict[str(ant) + '_freq'], adict[str(ant) + '_amp'], '.-')
            if A == 7:
                ax.set_xlabel('Freq')
                ax.set_ylabel('Sq Mag')

        pearsonr = scipy.stats.pearsonr(np.log10(adict['freq']), np.log10(adict['amp']))
        ant_corr.append(pearsonr[0])
        
        ax.set_title('Feed %i, Corr = %f' %(feed, pearsonr[0]))
 
    plt.draw()
    plt.savefig('ant_feed_corr/feeds_mult_ant.png', type='png')
        
    fig = plt.figure(2)
    plt.clf()
    plt.hist(ant_corr, bins=20)
    plt.text(-0.3, 1.5, 'Mean = %f' %np.average(ant_corr))
    plt.axis([-1, 0, 0, 2])
    plt.savefig('ant_feed_corr/feed_corr.png', type='png')
                

def beam_width_stats(wherecmd='', bins=50):
    """Plot histogram of beam width."""
    
    # Query data
    data = hextoolkit.querysql(['freq', 'x_width_el', 'x_width_az',
                                'y_width_el', 'y_width_az'], wherecmd=wherecmd)
    
    fghz = data['freq'] / 1000.
                
    xwe = 2 * data['x_width_el'] * fghz
    xwa = 2 * data['x_width_az'] * fghz
    ywe = 2 * data['y_width_el'] * fghz
    ywa = 2 * data['y_width_az'] * fghz
                
    fig = plt.figure(1, figsize=(12,12))
    plt.clf()
    
    print 'DATA: xwe, xwa, ywe, ywa'
    print 'MEDIAN:', np.median(xwe), np.median(xwa), np.median(ywe), np.median(ywa)
    print 'MEAN:  ', np.mean(xwe), np.mean(xwa), np.mean(ywe), np.mean(ywa)
    print ' ERR:  ', np.std(xwe) / np.sqrt(xwe.size), np.std(xwa) / np.sqrt(xwa.size), np.std(ywe) / np.sqrt(ywe.size), np.std(ywa) / np.sqrt(ywa.size)
    print 'STDEV: ', np.std(xwe), np.std(xwa), np.std(ywe), np.std(ywa)
    
    ax = fig.add_subplot(221)
    ax.hist(xwe, bins=bins)
    ax.set_title('X Width El')
    
    ax = fig.add_subplot(222)
    ax.hist(xwa, bins=bins)
    ax.set_title('X Width Az')
    
    ax = fig.add_subplot(223)
    ax.hist(ywe, bins=bins)
    ax.set_title('Y Width El')
    
    ax = fig.add_subplot(224)
    ax.hist(ywa, bins=bins)
    ax.set_title('Y Width Az')
 
 
def beam_width_vs_freq(bins=50):
    """Plot beam width by pol and direction."""
    
    # Query data
    data = hextoolkit.querysql(['freq', 'x_width_el', 'x_width_az',
                                'y_width_el', 'y_width_az'])
                
    freqlist = np.unique(data['freq'])
    xwe = []
    xwa = []
    ywe = []
    ywa = []
    xwe_uc = []
    xwa_uc = []
    ywe_uc = []
    ywa_uc = []
    
    for f in freqlist:
        fdata = data[np.where(data['freq'] == f)]
        
        fxwe = 2 * fdata['x_width_el'] * f / 1000.
        fxwa = 2 * fdata['x_width_az'] * f / 1000.
        fywe = 2 * fdata['y_width_el'] * f / 1000.
        fywa = 2 * fdata['y_width_az'] * f / 1000.
    
        xwe.append(np.mean(fxwe))
        xwa.append(np.mean(fxwa))
        ywe.append(np.mean(fywe))
        ywa.append(np.mean(fywa))
        xwe_uc.append(np.std(fxwe) / np.sqrt(fxwe.size))
        xwa_uc.append(np.std(fxwa) / np.sqrt(fxwa.size))
        ywe_uc.append(np.std(fywe) / np.sqrt(fywe.size))
        ywa_uc.append(np.std(fywa) / np.sqrt(fywa.size))
                
    fig = plt.figure(1, figsize=(14,12))
    plt.clf()
    
    for val,uc,name in [(xwe, xwe_uc, 'X Width El'), 
                        (xwa, xwa_uc, 'X Width Az'), 
                        (ywe, ywe_uc, 'Y Width El'), 
                        (ywa, ywa_uc, 'Y Width Az')]:
        plt.errorbar(freqlist / 1000., np.array(val), yerr=np.array(uc),
            fmt='.-', label=name)

    plt.legend(loc='lower right')
    plt.xlim([np.min(freqlist / 1000. - 1), np.max(freqlist / 1000. + 1)])
    plt.draw()
    
    plt.title('Beam Width vs. Frequency')
    plt.xlabel(r'$f_{GHz}$')
    plt.ylabel(r'$FWHM \times f_{GHz}$')
    
    plt.savefig('beamw_freq.png', type='png')




