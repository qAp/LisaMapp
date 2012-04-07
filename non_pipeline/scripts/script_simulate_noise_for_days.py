#!/usr/bin/env python
import os
import sys
import glob
import time

days = range( 1 , 10+1 )
stime = 8640.
Nb = 2
seed = 99 #This can be 'random', or a 'positive integer' ( as in a integer in a string.)
tsdir = 'data_tr1'


#Divide the days into batches so they can run in parallel
if days == None :
    print 'No days are selected.  Nothing to do.' ; sys.exit()
else :
    Ndays = len( days )
    if Ndays % Nb == 0 :
        Ndb = Ndays / Nb
        days_batches = [ days[ b*Ndb : (b+1)*Ndb  ]  for b in range( Nb ) ]
    elif Ndays % Nb > 0 :
        Ndb = Ndays / ( Nb-1 )
        Ndbl = Ndays % ( Nb-1 )
        days_batches = [ days[ b*Ndb : (b+1)*Ndb ] for b in range( Nb-1 ) ] + [ days[ (Nb-1)*Ndb : ] ]
    else :
        raise Exception , "Both the number of days and number of batches have to be postivie integer!"

#Create each batch's submit file and submit the jobs
for b in range( Nb ) :
    batch = b + 1
    print "Preparing and submitting batch %d" % batch
    if days_batches[ b ] == [] :
        print 'No days in this batch.  Nothing to do...' ; sys.exit()
    else :
        commands = [ ('./x_simulate_noise_for_days.py -d%d --stime %f --seed %s %s\n') %
                     ( day , stime , seed , tsdir ) for day in days_batches[ b ] ]
        submitname = 'x_simulate_noise_for_days_b%03d' % batch
        file = open( submitname + '.sub' , 'w' )
        file.writelines( [ '#!/bin/bash\n' ,
                           '#PBS -N %s\n' % ( submitname + '.sub' ) ,
                           '#PBS -o %s\n' % ( submitname + '.out' ) ,
                           '#PBS -j oe\n' ,
                           '#PBS -q compute\n' ,
                           '#PBS -l nodes=1:ppn=1\n' ,
                           '#PBS -l walltime=10:00:00\n' , 
                           'cd $PBS_O_WORKDIR\n' , '\n' ] +
                         commands + [ 'echo done' ] ) ; file.close()
        file = open( submitname + '.out' , 'w' ) ; file.write( 'dummy output' ) ; file.close()
        os.system( 'qsub %s' % ( submitname + '.sub' ) ) ; print 'done'
    time.sleep( 0.1 )



