#! /usr/bin/env python
## http://ast.noao.edu/sites/default/files/GMOS_Cookbook/Processing/PyrafProcIFU.html
## https://gmos-ifu-1-data-reduction-tutorial-gemini-iraf.readthedocs.io/en/latest/

# work_dir = '/lwk/xychen/AKARI_ULIRG/GMOS/observation/'
obs_dir = '/Users/eaarc/Documents/Gemini/observation/'

import sys
import copy
import numpy as np
from astropy.io import fits as fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as mplt
from sortedcontainers import SortedDict
from pyraf import iraf
from pyraf.iraf import gemini, gemtools, gmos, gfreduce, gfextract, gdisplay, rv, fitsutil
sys.path.append(obs_dir+'auxiliary/')
import fileSelect as fs
## task gdisplay_mod = "./gdisplay_mod.cl"
## force fl_erase = yes in iraf.gdisplay otherwise too slow

task lacos_spec_mod = "/Users/eaarc/Documents/Gemini/observation/auxiliary/lacos_spec_mod.cl"
# only use full path (not dir+'../'')

def nearest_file(obs, datalist) : 
    idx = np.abs(np.array([int(obs[-4:]) - int(x[-4:]) for x in datalist])).argmin()
    return datalist[idx]

# function of joining amplifers
def joinamp(input, exts, num=1) : 
    #iraf.copy( input+'.fits', input+'_join.fits' )
    iraf.delete( input+'_join.fits' )
    for ext in exts : # ['SCI', 'VAR', 'DQ'] : 
        for amp in ( [1] if (num==1) else [1, 5, 9] ) : 
            joinlist = ''
            for i in ( range(12) if (num==1) else range(4) ): 
                joinlist += input+'.fits['+ext+',' + str( amp+i ) + '],'
            iraf.delete( 'tmp_join.fits' )
            iraf.imjoin( joinlist[:-1], 'tmp_join.fits', '1', verbose='no' )
            #iraf.imcopy( 'tmp_join.fits', input+'_join.fits['+ext+','+str( amp )+',overwrite]')
            iraf.imcopy( 'tmp_join.fits', input+'_join.fits['+ext+','+str( amp )+']')

#######################
###### Related Files ######
#######################

## cd ./raw
## python2 obslog.py obsLog.sqlite3
dbFile = obs_dir+'raw/GN-2020B-Q-117/obsLog.sqlite3'
rawdir = obs_dir+'raw/GN-2020B-Q-117/'

qds_tmp = {'use_me':1,
       'CcdBin':'2 1',
       'DateObs':'FILL_LATER',
       'Instrument':'GMOS-N',
       'Disperser':'R150+_G5308',
       'AperMask':'IFU-2',
       'CentWave':-9999,
       'Object':'FILL_LATER',
       'RoI':'Full'
       }

qds_Sci = copy.deepcopy(qds_tmp)
qds_Sci['Object'] = 'J112657.76+163912.0'
qds_Sci_wide = copy.deepcopy(qds_Sci)
########
qds_Sci['DateObs'] = '2021-03-06:2021-03-06'
qds_Sci_wide['DateObs'] = '2021-03-05:2021-03-07'
Sci_comb = [ 'J1126_n1_comb' ]
########
# qds_Sci['DateObs'] = '2021-03-22:2021-03-22'
# qds_Sci_wide['DateObs'] = '2021-03-21:2021-03-23'
# Sci_comb = [ 'J1126_n2_comb' ]

###### Sci./Std. ######

Sci = {}
CentWaves = []
# for wave in [str(x) for x in range(500,900)] : 
for wave in ['680', '700', '720', '740', '760', '780'] : 
    qds_Sci['CentWave'] = int(wave)
    tmp = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qds_Sci), qds_Sci)
    if len( tmp ) >= 1 : 
        Sci[wave] = tmp
        CentWaves.append(wave)

Sci = SortedDict(Sci)
print ('Science observations:', Sci)

###### Bias ######

Bias = fs.fileListQuery(dbFile, fs.createQuery('bias', qds_Sci_wide), qds_Sci_wide)
print ('Bias observations:', Bias)
DateObs = qds_Sci['DateObs'].replace(":", "_").replace("-", "_")

###### Flat ######

Flat = {}
for wave in CentWaves:
    qds_Sci['CentWave'] = int(wave)
    files = fs.fileListQuery(dbFile, fs.createQuery('gcalFlat', qds_Sci), qds_Sci)
    tmp = []
    for obs in Sci[wave]:
        tmp.append(nearest_file(obs, files))
    Flat[wave] = sorted(list(set(tmp)))

Flat = SortedDict(Flat)
print ('Flat observations:', Flat)

for epoch in CentWaves:
    print(epoch)
    for obs in Sci[epoch]:
        print(obs, nearest_file(obs, Flat[epoch]))

###### CuAr ######

Arc = {}
for wave in CentWaves:
    qds_Sci_wide['CentWave'] = int(wave)
    files = fs.fileListQuery(dbFile, fs.createQuery('arc', qds_Sci_wide), qds_Sci_wide)
    tmp = []
    for obs in Sci[wave]:
        tmp.append(nearest_file(obs, files))
    Arc[wave] = sorted(list(set(tmp)))

Arc = SortedDict(Arc)
print ('CuAr observations:', Arc)

for epoch in CentWaves:
    print(epoch)
    for obs in Sci[epoch]:
        print(obs, nearest_file(obs, Arc[epoch]))

