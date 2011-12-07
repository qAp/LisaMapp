#!/usr/bin/env python
import os
import sys
import cPickle as cpkl


setup = {}

""" Where are the executables for simulation?"""
setup['execdir'] = ''#'/gpfs1/JC0311443/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'

""" TDI observables to simulate time-series for """
setup['tditype'] , setup['tdigen'] , setup['whichtdi'] = 'Michelson' , 'G2' , 'optimal'

""" Frequencies of the ORF multipole moments """
setup['f0'] , setup['df'] , setup['Nf'] = 0.0001 , 0.0001 , 4999


file = open( 'setup_get_expected_NSDs.pkl' , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()


