"""
Hex analysis toolkit.

Author: Keaton J. Burns <keaton.burns@gmail.com>
Affiliation: UC Berkeley

"""


import os
import os.path
import hextoolkit
import hexplot
import hexanalysis


plotdir = '/Users/kburns/Research/radio/plots/paper_images/'

def all():
    # Change to plot directory
    plotdir = '/Users/kburns/Research/radio/plots/paper_images/'
    os.chdir(plotdir)
    
    # Run Plots
    
    # Sumchisq distribution
    hexanalysis.abs_dist('sumchisq', 90, save=True)
        
    hexanalysis.beam_mag_angle()
    
    ha.magfreq_powerlaw_rev(save=True)




def dist_sumchisq():
    




def beam_shape():
    path = os.path.join(plotdir, 'beam_shape', 'scatter/')
    if not os.path.exists(path): os.makedirs(path)
    hexanalysis.beam_shape_plot('az_el', 'scatter', saveas=path+'az_el.png')
    hexanalysis.beam_shape_plot('mag_angle', 'scatter', saveas=path+'mag_angle.png')
    hexanalysis.beam_shape_plot('mag_ecc', 'scatter', saveas=path+'mag_ecc.png')    
    
    path = os.path.join(plotdir, 'beam_shape', 'density/')
    if not os.path.exists(path): os.makedirs(path)
    hexanalysis.beam_shape_plot('az_el', 'density', saveas=path+'az_el.png')
    hexanalysis.beam_shape_plot('mag_angle', 'density', saveas=path+'mag_angle.png')
    hexanalysis.beam_shape_plot('mag_ecc', 'density', saveas=path+'mac_ecc.png') 
    






