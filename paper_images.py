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


plotdir = '/Users/kburns/Research/radio/plots/'

def others()  
    ha.magfreq_powerlaw_rev(save=True)




def dist_sumchisq():
    path = os.path.join(plotdir, 'dist_sumchisq.png')
    hexanalysis.abs_dist('sumchisq', 90, saveas=path)

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
    






