#!/usr/bine/env python
import os
import sys
import cPickle as cpkl

"""
This sets up the simulation of data time-series lasting one day on different days of the year.
The option of preserving the seed (for random number generation) for repeated simulation is not enabled.
(It is not obvious how this can be done when simulating using Torque)
"""



setup = {}

""" Where are the executables for simulation?"""
setup['execdir'] = '' #'/gpfs1/JC0311443/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'

""" Path to the file containing the source distribution """
setup['Ppath'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/skymaps/library/quasi-onepixel/lon_263_lat_-35_plus_P00/P.pkl'

""" Directory of the ORF's multipole moments """
setup['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/myLISAmodule_stuff/EccentricInclined_eta0_0_xi0_0_sw_1_t0_0/tdiORF_SpHs/'

""" Name of the directory to save the simulated data in """
setup['tsdir'] = 'data/'

""" GW spectral slope """
setup['GWslope'] = 0

""" TDI observables to simulate time-series for """
setup['tditype'] , setup['tdigen'] , setup['whichtdi'] = 'Michelson' , 'G2' , 'optimal'

""" Maximum l for the multipole moments of ORF and Plm in the convolution """
setup['lmax'] = 20

""" Sample time of time-series """
setup['stime'] = 1.0

""" Frequencies of the ORF multipole moments """
setup['f0'] , setup['df'] , setup['Nf'] = 0.0001 , 0.0001 , 4999

""" On what days to simulate data? """
setup['days'] = range( 1 , 365+1 )

""" Number of batches to divided the days into """
setup['number of batches']  = 73

file = open( 'setup_signal_simulation.pkl' , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()


