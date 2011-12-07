#!/usr/bin/env python
import os
import sys
import glob

days = range( 1 , 2+1 )
stime = 1.
Nb = 1
seed = 'None'
tsdir = 'data_r1'



os.system( ('./x_simulate_noise_for_days.py ' + '-d%d '*len( days ) + '--stime %f --number_batches %d --seed %s %s')
           % tuple( days + [ stime , Nb , seed , tsdir ] ) ) == 0


