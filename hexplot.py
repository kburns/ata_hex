# hexplot.py
# Keaton Burns, University of California Berkeley, 09/28/2010
"""Tools for visualizing data from hex SQLite database"""


import numpy as np
import sqlite3
import matplotlib
matplotlib.use('pdf')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import hextoolkit


def hexplot(xdata, ydata, groupby=None, colorby=None, cmap=None, wherecmd='',
            saveas='squintplots.pdf', title=None, lines=False, errorbars=True, 
            exclude_flagged=True, xlim=None, ylim=None, return_data=False):
    """
    hexplot
    =======
    
    PURPOSE:
        Plots data from squint.db, grouped and colored as requested
    
    CALLING SEQUENCE:
        hexplot(xdata, ydata, groupby=None, colorby=None, wherecmd='',
                saveas='squintplots.pdf', title=None, lines=False, errorbars=True, 
                exclude_flagged=True, xlim=None, ylim=None)
    
    INPUTS:
        xdata       :=  tag for x-data in plots
        ydata       :=  tag for y-data in plots
        groupby     :=  group into separate plots by this tag
        colorby     :=  color by this tag
        cmap        :=  colormap for coloring dense data
        wherecmd    :=  'WHERE ...' command for specifying data in sql query
        saveas      :=  savename for pdf of plots
        title       :=  title to add to all plots
        lines       :=  whether to connect the plotted points with lines
                        (add ' ORDER BY ' to wherecmd to control line order)
        errorbars   :=  add errorbars when available
        exclude_flagged :=  remove datapoints which have been flagged
        xlim        :=  user-specified x-axis limits (xmin, xmax)
        ylim        :=  user-specified y-axis limits (ymin, ymax)
        return_data :=  return data array instead of plotting
    

    TAG LIST:
        In SQL database:
            'date'
            'source'
            'freq'
            'flux'
            'archsummpath'      (not physical)
            'rid'               (not physical)
            'antnum'
            'antname'
            'feed'
            'feedrev'               Decimal feed revision
            'sefd'
            'sumchisq'
            'squintaz'
            'squintaz_uc'
            'squintel'
            'squintel_uc'
            'x_width_az'
            'x_width_az_uc'
            'x_width_el'
            'x_width_el_uc'
            'y_width_az'
            'y_width_az_uc'
            'y_width_el'
            'y_width_el_uc'
            'flag'
            
        Derived from SQL database:
            'antfeed'
            'antfeedrev'
            'squintmag'
            'squintmag_uc'
            'squintangle'
            'squintangle_uc'
            
            
    TODO LIST:
        -look into interactive plotting
        -outlier identification (by uc?) & automatic flagging
        -more colors
        -mag vs feed by ant
           
    """
    

    # Query database
    data = hextoolkit.querysql([xdata, ydata, groupby, colorby], 
               wherecmd=wherecmd, exclude_flagged=exclude_flagged)
     
    # Setup pdf output
    pp = PdfPages(saveas)
    
    # Get number of and list of groups
    groupnum = 1
    if groupby != None:
        grouplist = np.unique(data[groupby])
        groupnum = np.size(grouplist)
        if groupnum > 100:
            print 'HEXPLOT: Requested grouping would yield over 100 plots'
            print 'HEXPLOT: Quitting...'
            pp.close()
            return
            
    # Get number of and list of coloring groups
    palette = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    colornum = 1
    if colorby != None:
        colorlist = np.unique(data[colorby])
        colornum = np.size(colorlist)
        if colornum > len(palette):
            cmap = True
            print 'HEXPLOT: Requested coloring would yield over 8 colors'
            print 'HEXPLOT: Using colormap instead'
            #pp.close()
            #return
                
    # Plot a separate figure for each group
    for i in xrange(groupnum):
        plt.figure(i+1, figsize=(9, 9))
        plt.clf()
        plt.subplots_adjust(left=0.125, right=0.9, bottom=0.125, top=0.9)
        
        # Pull out individual group, take everything if no grouping is specified
        if groupby != None: 
            igroup = np.where(data[groupby] == grouplist[i])
        else:
            igroup = np.where(data[xdata] == data[xdata])
            
        ixdata = data[xdata][igroup]
        iydata = data[ydata][igroup]

        # Determine errorbars as requested
        [xuc, yuc] = [None, None]
        
        if errorbars:
            if xdata in ['squintel', 'squintaz', 'squintmag', 'squintangle',
                         'x_width_el', 'x_width_az', 'y_width_el', 'y_width_az']:
                xuc = data[xdata + '_uc'][igroup]

            if ydata in ['squintel', 'squintaz', 'squintmag', 'squintangle',
                         'x_width_el', 'x_width_az', 'y_width_el', 'y_width_az']:
                yuc = data[ydata + '_uc'][igroup]

        # Plot the group
        if lines: linestyle = 'o-'
        else: linestyle = 'o'

        if colorby:
            if cmap:
                if errorbars:
                    print 'HEXPLOT: Cannot use colormap with errorbars'
                    print 'HEXPLOT: Errorbars will be supressed'
            
                if cmap == True:
                    cmap = matplotlib.cm.jet
                
                plt.scatter(ixdata, iydata, c=data[colorby][igroup], cmap=cmap)
    
            else:
                for j in xrange(colornum):
                    cgroup = np.where(data[colorby][igroup] == colorlist[j])
                    if np.size(cgroup) == 0: continue
                    
                    clabel = str(colorlist[j])                    
                    plotformat = palette[j] + linestyle
                    
                    [cxuc, cyuc] = [None, None]
                    if xuc != None: cxuc = xuc[cgroup]
                    if yuc != None: cyuc = yuc[cgroup]
                    
                    plt.errorbar(ixdata[cgroup], iydata[cgroup], xerr=cxuc, yerr=cyuc, 
                                 fmt=plotformat, label=clabel)
                             
        else:
            plt.errorbar(ixdata, iydata, xerr=xuc, yerr=yuc, fmt='b'+linestyle)
            
        # Add labels and lines at axes
        title_list = []
        if title != None: title_list.append(title)
        if groupby != None: title_list.append(groupby + ' = ' + str(grouplist[i]))
        plt.title(', '.join(title_list))
        plt.xlabel(xdata)
        plt.ylabel(ydata)
        plt.axhline(0, linestyle=':', color='k')
        plt.axvline(0, linestyle=':', color='k')
        if colorby: 
            if cmap:
                plt.colorbar()
            
            else:
                prop = matplotlib.font_manager.FontProperties(size='small')
                plt.legend(prop=prop, numpoints=1, title=colorby)
        
        # Data limits
        plotlimits = [np.min(ixdata), np.max(ixdata), np.min(iydata),np.max(iydata)]
        
        # Modify for particular tags:
        # Lower limit zero
        if ydata in ['squintaz_uc', 'squintel_uc', 'squintmag', 'squintmag_uc', 
                     'squintangle_uc', 'sefd', 'sumchisq']:
            plotlimits[2] = 0.0
        
        # Symmetric about zero
        if ydata in ['squintaz', 'squintel']:
            plotlimits[2] = -np.max(np.abs(plotlimits[2:4]))
            plotlimits[3] = np.max(np.abs(plotlimits[2:4]))
        if xdata in ['squintaz', 'squintel']:
            plotlimits[0] = -np.max(np.abs(plotlimits[0:2]))
            plotlimits[1] = np.max(np.abs(plotlimits[0:2]))
        
        # Specific values
        if ydata in ['squintangle']:
            plotlimits[2] = -np.pi
            plotlimits[3] = np.pi
        if xdata in ['squintangle']:
            plotlimits[0] = -np.pi
            plotlimits[1] = np.pi

        # Square up if comparing like quantities
        arcmin = ['squintaz', 'squintel', 'squintmag']
        arcmin_uc = ['squintaz_uc', 'squintel_uc', 'squintmag_uc']
        if xdata in arcmin and ydata in arcmin:
            plotlimits = [-np.max(plotlimits), np.max(plotlimits), 
                          -np.max(plotlimits), np.max(plotlimits)]
        elif xdata in arcmin_uc and ydata in arcmin_uc:
            plotlimits = [np.min(plotlimits), np.max(plotlimits), 
                          np.min(plotlimits), np.max(plotlimits)]
 
        # Allow user determination
        if xlim != None:
            plotlimits[0:2] = list(xlim)
        if ylim != None:
            plotlimits[2:] = list(ylim)
        
        # Pad plot limits on all sides
        pad = 0.05
        plotlimits[0] -= pad * (plotlimits[1] - plotlimits[0])
        plotlimits[1] += pad * (plotlimits[1] - plotlimits[0])
        plotlimits[2] -= pad * (plotlimits[3] - plotlimits[2])
        plotlimits[3] += pad * (plotlimits[3] - plotlimits[2])
        
        plt.axis(plotlimits)
            
        plt.draw()
        pp.savefig()
        
    pp.close()



