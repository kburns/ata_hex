# hexanalysis.py
# Keaton Burns, University of California Berkeley, 03/22/2011
"""
Data generators for plotting some stuff....
"""


import os
import numpy as np
import scipy
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
import hextoolkit


plotfolder = '~/Research/hex7/plots/'


def abs_dist(tag, percentile=95, flag=False, plot=True, save=False):
    """
    abs_dist
    ========
    
    PURPOSE:
        Determine cutoff percentile for specified tag, plot histogram, flag
    
    CALLING SEQUENCE:
        abs_dist(tag, percentile=95, flag=False, plot=True, save=False)
    
    INPUTS:
        tag         :=  tag for data to examine
        percentile  :=  percentile for cutoff, in range [0, 100]
        flag        :=  flagging number if flagging obs past cutoff is desired
        plot        :=  Boolean to create histogram plot
        save        :=  Boolean to save histogram plot

    TAG LIST:
            'squintaz'
            'squintaz_uc'
            'squintel'
            'squintel_uc'
            'sefd'
            'sumchisq'

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
        
        if save:
            plt.savefig('dist_' + tag + '.png', type='png')
    
    if flag:
        # Default flag value
        if flag is True:
            flag = 10
            
        hextoolkit.flagsql('WHERE oid IN (SELECT oid FROM runs natural join obs WHERE abs(' + tag + ') >= %f)' %cutoff, num=flag)
     
 
def flagging(analyze=False, saveas=False):
    """Store flag settings"""
    
    # Reset all flags
    hextoolkit.flagsql(num=0, add=False)
    
    # Use only 3c48
    #hextoolkit.flagsql(wherecmd="WHERE rid IN (SELECT rid FROM runs WHERE source != '3c48')", num=2**14)

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
                            
    #     print np.sum(flaggrid, 0)
    #                            
    #     print [row[0] for row in flaglist]
    #     for n in xrange(N):
    #         print Alist[n], flaggrid[n]
            
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
            if saveas is True: saveas = 'current_flagging'
            plt.savefig(saveas + '.png')
        
    
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
    
def magfreq_powerlaw_fit(mincount=1, bins = 25, fit_of_medians=False, rev=False):
    """Compute power laws of squint vs mag dependence for each antfeedrev"""

    # Query data
    data = hextoolkit.querysql(['antfeedrev','date','squintmag','freq'])
    
    antfeedrevs = np.unique(data['antfeedrev'])
    freqlist = np.unique(data['freq'])
    daylist = np.unique(np.round(data['date']))

    slopes = []
    powers = []

    # Examine each antfeedred
    for afr in antfeedrevs:
        if rev:
            if np.round(np.mod(afr, 1), 1) != rev: continue

        idata = data[np.where(data['antfeedrev'] == afr)]
        
        # Take best fit of the median squintmag for each frequency
        if fit_of_medians:
            ifreqs = []
            imedians = []
            
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
                
                ifreqs.append(f)
                imedians.append(ifmag)
                
            if len(ifreqs) < 2: continue
            
            # MIN AMPLITUDE:
            ave_amp = np.average(imedians)
            if ave_amp > 4:
                print afr, ave_amp
            else:
                continue
                
            # Compute power law and linear slopes    
            logfit = np.polyfit(np.log10(ifreqs), np.log10(imedians), 1)
            powers.append(logfit[0])
            
            linfit = np.polyfit(ifreqs, imedians, 1)
            slopes.append(linfit[0])
        
        # Take median of the best fit of each day
        else:
            logfits = []
            linfits = []
            
            # Examine each day
            for d in daylist:
                iddata = idata[np.where(np.round(idata['date']) == d)]
                # Skip if no points for this antfeedrev in this run
                if iddata.size < 2: 
                    continue
                
                # Compute power law and linear slopes    
                logfit = np.polyfit(np.log10(iddata['freq']), 
                                    np.log10(iddata['squintmag']), 1)
                logfits.append(logfit[0])
                
                linfit = np.polyfit(iddata['freq'], iddata['squintmag'], 1)
                linfits.append(linfit[0])
                
            if len(logfits) < mincount: 
                continue
            
            # Get statistics of data closest to median
            imedianlog = np.median(logfits)
            logmedindex = (np.abs(logfits - imedianlog)).argmin()
            
            imedianlin = np.median(linfits)
            linmedindex = (np.abs(linfits - imedianlin)).argmin()
            
            #powers.append(logfits[logmedindex])
            #slopes.append(linfits[linmedindex])
            for p in logfits:
                powers.append(p)
            
    return np.array(powers), np.array(slopes)
    
        
def magfreq_powerlaw_rev(mincount=1, bins=50, save=False, fit_of_medians=False):
    """Compute power laws by rev of squint vs mag dependence for each antfeedrev
    
    -look at large amplitudes
    -include uncertainties in median
    -test hypothesis -1 fit
    -goodness of fit
    
    """

    rev = [0.1, 0.2, 0.3]

    # Plot histograms of fits
    fig = plt.figure(1, figsize=(10,10))
    fig.clf()

    for i in xrange(len(rev)):
        powers, slopes = magfreq_powerlaw_fit(mincount=mincount, fit_of_medians=fit_of_medians, rev=rev[i])
        
        # Disregard outliers
        outliers = (np.abs(powers) >= 5)
        print '%i outliers of %i' %(outliers.sum(), powers.size)
        powers = powers[~outliers]
                        
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
    
def beam_width_scatter():

    # Query data
    data = hextoolkit.querysql(['freq', 'x_width_el', 'x_width_az',
                                'y_width_el', 'y_width_az'])
                                
    fghz = data['freq'] / 1000.
    xwe = data['x_width_el']
    xwe_uc = data['x_width_el_uc']
    xwa = data['x_width_az']
    xwa_uc = data['x_width_az_uc']
    ywe = data['y_width_el']
    ywe_uc = data['y_width_el_uc']
    ywa = data['y_width_az']
    ywa_uc = data['y_width_az_uc']
    
    # Definitions for the axes
    nullfmt = matplotlib.ticker.NullFormatter()
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    for pol in ['x', 'y']:

        plt.figure(1, figsize=(8,8))
        plt.clf()
        
        if pol == 'x':
            x = 2 * xwa * fghz
            xuc = 2 * xwa_uc * fghz
            y = 2 * xwe * fghz
            yuc = 2 * xwe_uc * fghz
        else:
            x = 2 * ywa * fghz
            xuc = 2 * ywa_uc * fghz
            y = 2 * ywe * fghz
            yuc = 2 * ywe_uc * fghz
    
        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        
        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        # the scatter plot:
        axScatter.errorbar(x, y, xerr=xuc, yerr=yuc, fmt='o')
        axScatter.set_xlabel(r'$\rm{AZ FWHM}_{%s} \times f_{GHz}$' %pol)
        axScatter.set_ylabel(r'$\rm{EL FWHM}_{%s} \times f_{GHz}$' %pol)
        
        # now determine nice limits by hand:
        binwidth = 0.1
        xymax = np.max([np.max(x), np.max(y)])
        lim = (int(xymax/binwidth) + 1) * binwidth
        
        axScatter.set_xlim((0, lim))
        axScatter.set_ylim((0, lim))
        
        bins = np.arange(0, lim + binwidth, binwidth)
        axHistx.hist(x, bins=bins, log=True)
        axHisty.hist(y, bins=bins, orientation='horizontal', log=True)
        
        axHistx.set_xlim(axScatter.get_xlim())
        axHisty.set_ylim(axScatter.get_ylim())

        plt.savefig(pol + '_widths.png', type='png')
        
        
        
def beam_mag_angle():
    
    # Query data
    data = hextoolkit.querysql(['freq', 'xmag', 'xangle', 'ymag', 'yangle'])
                                
    fghz = data['freq'] / 1000.
    xm = data['xmag'] * 2 * fghz
    xm_uc = data['xmag_uc'] * 2 * fghz
    xa = data['xangle'] * 180. / np.pi
    xa_uc = data['xangle_uc'] * 180. / np.pi
    ym = data['ymag'] * 2 * fghz
    ym_uc = data['ymag_uc'] * 2 * fghz
    ya = data['yangle'] * 180. / np.pi
    ya_uc = data['yangle_uc'] * 180. / np.pi
    
    print 'DATA: x mag, y mag, x angle, y angle'
    print 'MEDIAN:', np.median(xm), np.median(ym), np.median(xa), np.median(ya)
    print 'MEAN:  ', np.mean(xm), np.mean(ym), np.mean(xa), np.mean(ya)
    print ' ERR:  ', np.std(xm) / np.sqrt(xm.size), np.std(ym) / np.sqrt(ym.size), np.std(xa) / np.sqrt(xa.size), np.std(ya) / np.sqrt(ya.size)
    print 'STDEV: ', np.std(xm), np.std(ym), np.std(xa), np.std(ya)
    
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
            x, xuc, y, yuc = xm, xm_uc, xa, xa_uc
        else:
            x, xuc, y, yuc = ym, ym_uc, ya, ya_uc
                
        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        
        # Remove histogram ticklabels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        # Scatter plot
        axScatter.errorbar(x, y, xerr=xuc, yerr=yuc, fmt='o')
        axScatter.axhline(45, c='r')
        axScatter.axvline(3.5, c='r')
        axScatter.set_xlabel(r'$\rm{%s\,Beam\,Size} \, \times \, f_{GHz}$' %pol)
        axScatter.set_ylabel(r'$\rm{%s\,Beam\,Angle\;(degrees)}$' %pol)
        
        # Histograms
        bins = 40
        xlim = axScatter.get_xlim()
        ylim = axScatter.get_ylim()
        axHistx.hist(x, bins=bins, log=True)
        axHistx.axvline(3.5, c='r')
        axHisty.hist(y, bins=bins, orientation='horizontal', log=True)
        axHisty.axhline(45, c='r')
        
        axHistx.set_xlim(xlim)
        axHisty.set_ylim(ylim)

        plt.savefig(pol + '_magangle.png', type='png')
        
        
def beam_mag_ecc():
    
    # Query data
    data = hextoolkit.querysql(['freq', 'xmag', 'ymag'])
                                
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
    
    # Definitions for the axes
    nullfmt = matplotlib.ticker.NullFormatter()
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    for pol in ['x', 'y']:

        plt.figure(1, figsize=(8,8))
        plt.clf()
        
        if pol == 'x':
            x = 2 * xm * fghz
            xuc = 2 * xm_uc * fghz
            y = xecc
            yuc = xecc_uc
        else:
            x = 2 * ym * fghz
            xuc = 2 * ym_uc * fghz
            y = yecc
            yuc = yecc_uc
    
        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        
        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        # the scatter plot:
        axScatter.errorbar(x, y, xerr=xuc, yerr=yuc, fmt='o')
        axScatter.set_xlabel(r'$\rm{%s\,Beam\,Size} \, \times \, f_{GHz}$' %pol)
        axScatter.set_ylabel(r'$\rm{%s\,Eccentricity}$' %pol)
        axScatter.set_ylim([-1, 1])
        print np.max(xuc)
        print np.max(yuc)
        
        # now determine nice limits by hand:
        bins = 40
        axHistx.hist(x, bins=bins, log=True)
        axHisty.hist(y, bins=bins, orientation='horizontal', log=True)
        
        axHistx.set_xlim(axScatter.get_xlim())
        axHisty.set_ylim(axScatter.get_ylim())

        plt.savefig(pol + '_magecc.png', type='png')
        
        
def squintmag_time_evolution(rev=False):
    
    # Query data
    data = hextoolkit.querysql(['antfeedrev','date','squintmag','freq'])
    
    antfeedrevs = np.unique(data['antfeedrev'])
    freqlist = np.unique(data['freq'])
    daylist = np.unique(np.round(data['date']))

    deviations = []
    
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
            
            # Compute standard deviation as a fraction of the mean
            squint = ifdata['squintmag']
            squint /= squint.mean()
            
            deviations.append(np.std(squint))

    return np.array(deviations)
    
def squintmag_time_rev(bins=50, save=False):

    rev = [0.1, 0.2, 0.3]

    # Plot histograms of fits
    fig = plt.figure(1, figsize=(10,10))
    fig.clf()

    for i in xrange(len(rev)):
        dev = squintmag_time_evolution(rev=rev[i])
                        
        # Power law histogram
        ax = fig.add_subplot(3, 1, i + 1)
        ax.hist(dev, bins=np.linspace(0, 2, bins))
        ax.set_ylabel(r'$\rm{Count}$')
        ax.text(0.75, 0.8, r'$\rm{Feed Rev} = %.1f$' %rev[i], transform=ax.transAxes)
        ax.text(0.75, 0.7, r'$\rm{Mean} = %.2f \pm %.2f$' %(np.mean(dev), np.std(dev) / np.sqrt(dev.size)), transform=ax.transAxes)
        #ax.text(0.05, 0.7, r'$\rm{StDev} = %.2f$' %np.std(dev), transform=ax.transAxes)

        if i == 2:
            ax.set_xlabel(r'$\rm{Squint \ Mag. \ Fractional \ Standard \ Devation}$')
            
    if save:
        if save == True: save = 'squinttime_rev.png'
        plt.savefig(save, type='png')
        

