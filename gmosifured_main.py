
# Must keep a blank line after for loop with a single command line.

###############################################################
##############################################################

###### Bias process ######

gemtools.gemextn.unlearn()    # Disarm a bug in gbias
gmos.gbias.unlearn()

procBias = 'mcBias_'+DateObs
iraf.imdelete(procBias)
for obs in Bias : iraf.imdelete('g'+obs)
# The str.join() function is needed to transform a python list into
# a string of comman-separated files that IRAF can understand.
gmos.gbias( ','.join([str(obs) for obs in Bias]), procBias, rawpath = rawdir, fl_vardq = 'yes', logfile = 'gbias.log' )

##############################################################
##############################################################

###### Verify the MDF ######

for wave in CentWaves:
    for obs in Flat[wave] :
        print(obs)
        iraf.copy('gmos$data/gnifu_slits_mdf.fits', 'MDF_'+obs+'.fits', verbose='no')

# prepare reduced Flat and Sci. data and attach temp MDF file
for wave in CentWaves :
    for obs in Flat[wave] :
        print(obs)
        iraf.imdelete('g'+obs)
        iraf.imdelete('rg'+obs)
        iraf.gfreduce(obs, \
            rawpath=rawdir, \
            mdffile='MDF_'+nearest_file(obs,Flat[wave])+'.fits', mdfdir='./', \
            bias='mcBias_'+DateObs, \
            fl_extract='no', fl_over='yes', fl_trim='yes', \
            fl_fluxcal='no', fl_gscrrej='no', fl_wavtran='no', fl_skysub='no', \
            fl_inter='no', fl_vardq='yes', slits='both', logfile = 'gfreduce.log')
            #l_addmdf='no', fl_bias='no', \

# require all the MDF to be the same for following combine
for wave in CentWaves :
    for obs in Flat[wave] :
        print(obs)
        iraf.tcalc('rg'+obs+'.fits[MDF]', 'BEAM', 'if NO == 750 then -1 else BEAM')
        iraf.tcalc('rg'+obs+'.fits[MDF]', 'BEAM', 'if NO == 795 then -1 else BEAM')
# tread('rg'+Flat_specified['680']+'.fits[MDF]') # view MDF, ^D quit to quit
# tstat('rg'+Flat_specified['680']+'.fits[MDF]', 'beam', rows='1-750', nrows=1) 
# tintegrate('rg'+Flat_specified['680']+'.fits[MDF]', 'beam', '')
# # 1, 775, 1152, 1155, already marked
# # 795 is marginnally worked
# for slit 1, the number of aperture is 748
# for slit 2, the number of aperture is 746 (if fiber 795 is wrong) or 747
## repeat the above MDF-edting and gfextract for checking

##############################################################
##############################################################

###### Extract Flat for tracing fibers of Arc/Sci. ######

prefix = 'rg'
for wave in CentWaves :
    for obs in Flat[wave] :
        print(obs)
        os.system( 'rm -v database/ap'+'e'+prefix+obs+'*' )
        iraf.imdelete( 'e'+prefix+obs )
        iraf.gfextract( prefix+obs, nsum=2000, recenter='yes', trace='yes', weights='none', fl_vardq='yes')

# check if the numbers match for both slits
for wave in CentWaves :
    for obs in Flat[wave] :
        print(wave, obs)
        tmp = fits.open( 'erg'+obs+'.fits' )
        for slit in [1,2] : 
            apnum = tmp['sci',slit].shape[0]
            column = [ float(x) for x in tmp['sci',slit].header[ 'DETSEC' ][1:-1].split(',')[0].split(':') ]
            col_apall_in = (column[1]-column[0]+1) / float( tmp['sci',slit].header[ 'CCDSUM' ].split()[0] )
            col_apall_out = float( tmp['sci',slit].header[ 'NAXIS1' ] )
            goodApall = col_apall_in == col_apall_out
            apIDs = np.array( ['APNUM'+str(x) for x in range(1, apnum+1)] )
            line0 = np.array( [ float(tmp['sci',slit].header[ apID ].split(' ')[2]) for apID in apIDs ] )
            line1 = np.array( [ float(tmp['sci',slit].header[ apID ].split(' ')[3]) for apID in apIDs ] )
            apSteps = line0[1:] - line1[:-1]
            isNormal = ( apSteps ) >= -4.5/2
            isPositive = ( apSteps ) >= 0
            #print( apIDs[1:][~isPositive], apSteps[~isPositive] )
            if sum(isNormal)+1 == apnum :
                print( slit, 'Normal', sum(isPositive)+1, goodApall, col_apall_in, col_apall_out )
            else : 
                print( slit, 'Error', sum(isPositive)+1, goodApall, col_apall_in, col_apall_out )
                print( 'APNUM1', tmp['sci',slit].header[ 'APNUM1' ] )
                for (apID0, apID1, apStep) in zip(apIDs[:-1][~isPositive], apIDs[1:][~isPositive], apSteps[~isPositive] ) : 
                    print( apID0, tmp['sci',slit].header[ apID0 ], apStep )
                    print( apID1, tmp['sci',slit].header[ apID1 ], apStep ) 
# common issue, automatically skipped:
# slit 1: 251/252, 643/644/645, 696/697, 741/742, apSteps ~ [-0.5,0]
# slit 2: 73/74, 825/826, 1480/1481, apSteps ~ [-0.5,0]

# correct with either method 1 or 2
# method 1, use the other good one as reference.
# the reference should have similar or larger apall columns
###################################################
#### Depends on different data ####
reference, files = 'N20210306S0054', [ 'N20210306S0050' ] # n1
reference, files = 'N20210322S0091', [ 'N20210322S0087' ] # n2
###################################################

for obs in files :
    print(obs)
    os.system( 'rm -v database/ap'+'e'+prefix+obs+'*' )
    iraf.imdelete( 'e'+prefix+obs )
    iraf.gfextract( prefix+obs, reference='e'+prefix+reference, recenter='yes', trace='yes', weights='none', fl_vardq='yes')
    # still keep recenter='yes', trace='yes' to allow some variance

########

# method 2, if no good reference, mannually create the reference
# check and repeat extract interactively
# the identification of fibres relies on nsum (i.e., which column(s))
# keep the loop format otherwise some input errors
for wave in CentWaves[0:1] :
    for obs in Flat[wave][0:1] :
        print(obs)
        os.system('rm -v database/aperg'+obs+'*')
        iraf.imdelete('erg'+obs)
        iraf.gfextract('rg'+obs, exslits='red', nsum=100, recenter='yes', trace='yes', weights='none', fl_inter='yes')
        os.system('rm -v database/aperg'+obs+'*')
        iraf.imdelete('erg'+obs+'*')
        iraf.gfextract('rg'+obs, exslits='blue', nsum=100, recenter='yes', trace='yes', weights='none', fl_inter='yes')
        os.system('rm -v database/aperg'+obs+'*')
        iraf.imdelete('erg'+obs+'*')
## - Answer 'yes' to the "Find apertures" question.
## - Answer 'yes' to the "Edit apertures" question.
## Each bundle is 50 fibers. 
## To zoom in, the interactive commands are “w”, 
## followed by “e-e” to defined the lower-left and upper-right corners. 
## “w” and “a” to zoom back out.
## - Type "q" to get out of the plot
## - Answer "NO" (uppercase) to all the questions.

# # generate the tracing file in database for two slits
# ###################################################
# #### Depends on different data ####
# nsums = [2000,100] # n1
# nsums = [100,100] # n2
# ###################################################
# obs = Flat[CentWaves[0]][0]
# os.system( 'rm -v database/ap'+'e'+prefix+obs+'*' )
# iraf.imdelete( 'e'+prefix+obs )
# iraf.gfextract( prefix+obs, exslits='red', nsum=nsums[0], recenter='yes', trace='yes', weights='none')
# iraf.imdelete('e'+prefix+obs+'*')
# iraf.gfextract( prefix+obs, exslits='blue', nsum=nsums[1], recenter='yes', trace='yes', weights='none')
# iraf.imdelete('e'+prefix+obs+'*')
# # trace by self
# iraf.gfextract( prefix+obs, reference='e'+prefix+obs, recenter='yes', trace='yes', weights='none', fl_vardq='yes')
# # go back to method 1 for other flats

##############################################################
##############################################################

###### Preparation with corrected MDF from Flat ######

# update MDF files for usage of Arc/Sci/Std data
for wave in CentWaves :
    for obs in Flat[wave] :
        print(obs)
        iraf.delete('MDF_'+obs+'.fits')
        iraf.tcopy('rg'+obs+'.fits[MDF]', 'MDF_'+obs+'.fits', verbose='yes')
# tread('MDF_'+Flat_specified['740']+'.fits') # view MDF, ^D quit to quit
# OR directly modify in Arc/Sci/Std file:
# iraf.tcalc('rg'+Std['680'][0]+'.fits[MDF]', 'BEAM', 'if NO == 795 then -1 else BEAM')

for wave in CentWaves :
    for obs in Arc[wave] + Sci[wave] :
        print(obs)
        iraf.imdelete('g'+obs)
        iraf.imdelete('rg'+obs)
        iraf.gfreduce(obs, \
            rawpath=rawdir, \
            mdffile='MDF_'+nearest_file(obs,Flat[wave])+'.fits', mdfdir='./', \
            bias='mcBias_'+DateObs, \
            fl_extract='no', fl_over='yes', fl_trim='yes', \
            fl_fluxcal='no', fl_gscrrej='no', fl_wavtran='no', fl_skysub='no', \
            fl_inter='no', fl_vardq='yes', slits='both', logfile = 'gfreduce.log')
            #l_addmdf='no', fl_bias='no', \

###### Extract Arc and Sci tracing by ergFlat ######

prefix = 'rg'
for wave in CentWaves :
    for obs in Arc[wave] + Sci[wave] :
        print(obs)
        os.system( 'rm -v database/ap'+'e'+prefix+obs+'*' )
        iraf.imdelete( 'e'+prefix+obs )
        iraf.gfextract( prefix+obs, reference='erg'+nearest_file(obs,Flat[wave]), recenter='no', trace='no', weights='none', fl_vardq='yes')

# check if image column matches DETSEC 
for wave in CentWaves :
    for obs in Arc[wave] + Sci[wave] : 
        print(wave, obs)
        tmp = fits.open( 'erg'+obs+'.fits' )
        for (slit, apnum) in zip([1,2], [748,746]) : 
            column = [ float(x) for x in tmp['sci',slit].header[ 'DETSEC' ][1:-1].split(',')[0].split(':') ]
            col_apall_in = (column[1]-column[0]+1) / float( tmp['sci',slit].header[ 'CCDSUM' ].split()[0] )
            col_apall_out = float( tmp['sci',slit].header[ 'NAXIS1' ] )
            goodApall = col_apall_in == col_apall_out
            print( 'spec_col:', col_apall_in, col_apall_out, goodApall )  

# # check the cross-matching between fiber ID and its locaiton
# for wave in CentWaves :
#     print(wave)
#     for obs in Sci[wave] :
#         print(obs)
#         iraf.gfdisplay( 'erg'+obs, 1, z1=15e4, z2=60e4, ztrans='linear' )

##############################################################
##############################################################

###### Wavelength Solution ######

prefix = 'rg'
for wave in CentWaves :
    print(wave)
    for obs in Arc[wave] :
        print(obs)
        os.system( 'rm -v database/id'+'e'+prefix+obs+'*' )
        iraf.gswavelength( 'e'+prefix+obs, \
            fl_inter='yes', nsum=200, nlost=10, ntarget=15, threshold=25, \
            coordlis=obs_dir+'auxiliary/wavecalib/GCALcuar_select.dat', logfile = 'gfreduce.log')
#         coordlis='gmos$data/GCALcuar.dat', logfile = 'gfreduce.log')
## - Type "l" to automatically fit the lines the tool can fit now.
## - Type "d" to remove weak identifications, only keep 
##   1) 4764.862 / 4879.86 (doublets, backup); 
##   2) 6965.43 (init) / 7067.32 (init)
##   3) 7272.935 / 7635.106 / 7723.98 
##   4) 7948.176 / 8264.522 / 8521.443 / 8667.944
##   5) 9122.967 / 9224.499 / 9657.78 / 9784.50 
##   6) 10470.0535
## - To mark a line, "m" with the cursor on the line.
## - Enter the wavelength in the text box if INDEF and press "Return".
## - Type "f" to get a fit with the list of lines.
## - Quit the "fit window" with "q".  (Just once!)
## - Type "q" to return to the "spectrum window".
## - Answer "NO" (uppercase) to "Fit dispersion function interactively".
# typical rms ~< 0.3
# CentWave 0th 
# 680 6469
# 700 6463
# 720 6457
# 740 6451
# 760 6446
# 780 6440
		
##############################################################
##############################################################

###### Find the gaps using Flat ######

# used for CR rejection by block, and scattered light subtraction
prefix = 'rg'
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] :
        print(obs)
        iraf.delete('blkmask_'+prefix+obs)
        iraf.gffindblocks( prefix+obs, 'erg'+obs, 'blkmask_'+prefix+obs)

# check the differences between the gap files
blkmask_ref = np.genfromtxt( 'blkmask_'+prefix+Flat[CentWaves[0]][0], skip_header = 0, names = [ 'col0', 'col1', 'line0', 'line1' ], dtype = 'i8, i8, i8, i8' ) 
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] :
        print(obs)
        blkmask = np.genfromtxt( 'blkmask_'+prefix+obs, skip_header = 0, names = [ 'col0', 'col1', 'line0', 'line1' ], dtype = 'i8, i8, i8, i8' ) 
        print('line0 diff: ', blkmask['line0'] - blkmask_ref['line0'])
        print('line1 diff: ', blkmask['line1'] - blkmask_ref['line1'])
# < +/-1 tolerance? 

# edit blkmask if required
# iraf.copy('blkmask_prgN20210516S0051_comb', 'blkmask_prgN20210516S0047_comb')
# blkmask = np.genfromtxt( 'blkmask_'+prefix+Flat_specified[epoch], skip_header = 0, names = [ 'col0', 'col1', 'line0', 'line1' ], dtype = 'i8, i8, i8, i8' ) 
# blkmask['line0'][0], blkmask['line1'][0] = 73, 77 # width around 4--5
# np.savetxt( blkmask_file, blkmask, fmt='%1i %3i %4i %4i')

##############################################################
##############################################################

###### CR rejection for each fiber region of Sci. (skip Flat) ######

prefix = 'rg'

for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] :
        print(obs)
        iraf.delete( 'x'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'x'+prefix+obs+'.fits' )
        ########################################
        # read gain and read noise from header
        tmp = fits.open( prefix+obs+'.fits' )
        isSciObs = tmp[0].header['object'] == qds_Sci['Object']
        RDNOISE = np.array( [ tmp[amp].header['RDNOISE'] for amp in range(1, 12+1) ] )
        GAIN = np.array( [ tmp[amp].header['GAIN'] for amp in range(1, 12+1) ] )
        # calculate mean (used below) and std, and check the ratio of std/mean
        print('GAIN mean, std, std/mean', np.mean(GAIN), np.std(GAIN), np.std(GAIN)/np.mean(GAIN))
        print('RDNOISE mean, std, std/mean', np.mean(RDNOISE), np.std(RDNOISE), np.std(RDNOISE)/np.mean(RDNOISE))
        # use the block mask from Flat 
        blkmask = np.genfromtxt( 'blkmask_rg'+nearest_file(obs,Flat[wave]), skip_header = 0, names = [ 'col0', 'col1', 'line0', 'line1' ], dtype = 'i8, i8, i8, i8' ) 
        block, line0, line1 = [], [], []
        # record the range of each gap
        for ind in range(0, len(blkmask)):
            block.append(0)
            line0.append(blkmask['line0'][ind])
            line1.append(blkmask['line1'][ind])
        # record the range of each block
        for ind in range(0, len(blkmask)-1):
            block.append(1)
            line0.append(blkmask['line1'][ind]+1)
            line1.append(blkmask['line0'][ind+1]-1)
        ########################################
        for ind in range(0, len(block)):
            print(line0[ind], line1[ind])
            joinlist_SCI, joinlist_DQ = '', ''
            for amp in range(1, 12+1) : 
                ypos = str(line0[ind])+':'+str(line1[ind])
                joinlist_SCI += prefix+obs+'.fits[SCI,'+str(amp)+'][*,'+ypos+'],'
                joinlist_DQ += prefix+obs+'.fits[DQ,'+str(amp)+'][*,'+ypos+'],'
            # create amp-joined image in the block
            iraf.imdelete( 'tmpblock_joinamp_SCI' )
            iraf.imjoin( joinlist_SCI[:-1], 'tmpblock_joinamp_SCI.fits', '1', verbose='no' )
            iraf.imdelete( 'tmpblock_joinamp_DQ' )
            iraf.imjoin( joinlist_DQ[:-1], 'tmpblock_joinamp_DQ.fits', '1', verbose='no' )
            ########################################
            # use modified lacos_spec to obtain clean SCI and mask
            iraf.imdelete( 'tmpclean1' )
            iraf.imdelete( 'tmpmask1' )
            iraf.lacos_spec_mod('tmpblock_joinamp_SCI', 'tmpclean1', 'tmpmask1', gain=np.mean(GAIN), readn=np.mean(RDNOISE), xorder=20*block[ind], yorder=3*block[ind], sigclip=10, sigfrac=0.5, objlim=0.5, niter=4, verbose='no', keeptmp='no') # xorder=(20 if isSciObs else 120)*block[ind]
            # use objlim=0.5 to remove very bright CR, which shows high obj.flux (high-underlying)
            # use sigclip=10 to protect bright sky lines (default sigclip=4.5); CR has higher sigma
            # relatively faint CR (but with sharp edge) remains in this step
            # only fit for block (block=1); no fit for gap (block=0)
            ########################################
            iraf.imdelete( 'tmpclean2' )
            iraf.imdelete( 'tmpmask2' )
            iraf.lacos_spec_mod('tmpclean1', 'tmpclean2', 'tmpmask2', gain=np.mean(GAIN), readn=np.mean(RDNOISE), xorder=5*block[ind], yorder=3*block[ind], sigclip=4.5, sigfrac=0.5, objlim=5, niter=4, verbose='no', keeptmp='no') # , keepmod='yes', outmod='tmpfitmod2'
            # use lower xorder to avoid fitting CR in x-axis (usually at the edge)
            # use sigclip=4.5 (default) to identify relatively faint CR
            # use objlim=5 to protect bright sky lines (default objlim=2)
            # only fit for block (block=1); no fit for gap (block=0)
            ########################################
            # caiculate the CR flux and recored into VAR
            iraf.imdelete( 'tmpvar2' )
            iraf.imarith( 'tmpblock_joinamp_SCI', '-', 'tmpclean2', 'tmpvar2' )
            # add mask of two step
            iraf.imarith( 'tmpmask1', '+', 'tmpmask2', 'tmpmask2' )
            iraf.imreplace( 'tmpmask2', 1, lower=1 )
            # scale tmpmask2 to DQ=8 and add into original DQ 
            iraf.imarith( 'tmpmask2', '*', 8, 'tmpmask2' )
            iraf.imarith( 'tmpblock_joinamp_DQ', '+', 'tmpmask2', 'tmpmask2' )
            ########################################
            # copy the corrected block back to original MEF format file
            for amp in range(1, 12+1) : 
                xpos = str(1+256*(amp-1))+':'+str(amp*256)
                ypos = str(line0[ind])+':'+str(line1[ind])
                iraf.imcopy( 'tmpclean2.fits['+xpos+',*]', 'x'+prefix+obs+'.fits[SCI,'+str(amp)+'][*,'+ypos+']' )
                iraf.imcopy( 'tmpvar2.fits['+xpos+',*]', 'x'+prefix+obs+'.fits[VAR,'+str(amp)+'][*,'+ypos+']', verbose='no' )
                iraf.imcopy( 'tmpmask2.fits['+xpos+',*]', 'x'+prefix+obs+'.fits[DQ,'+str(amp)+'][*,'+ypos+']', verbose='no' )

# delete tmp files
os.system( ' rm -v tmp* ' )

# just copy file for Flat
prefix = 'rg'
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave]:
        print(obs)
        iraf.delete( 'x'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'x'+prefix+obs+'.fits' )

#############
# check results:
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] :
        print(obs)
        iraf.delete( 'x'+prefix+obs+'_join.fits' )
        tmp = joinamp( 'x'+prefix+obs, ['SCI', 'VAR', 'DQ'] )
        iraf.display( 'x'+prefix+obs+'_join.fits[SCI]', 1 )
        iraf.display( 'x'+prefix+obs+'_join.fits[VAR]', 2, zscale='no', zrange='no', z1=-5, z2=100)
        raw_input( ) # pause

##############################################################
##############################################################

###### Bad pixel mask (Flat and Sci.) and fix for (only Flat) ######

# Outliers (edges, 0th-spec, CR) in Flat will be removed in this step. 
##############
# DQ convention: 
# 0: Good pixel
# 1: (Not found) Detector defect causing a bad pixel (hot or dead)
# 2: (Not found) Non-linear regime (may not be used for all instruments)
# 4: Saturated --> copy/grow to DQ=64
# 8: Cosmic-ray-hit
# 16: (Not found) No data
# 32: (Added) CCD edges (only for Flat)
# 64: (Added) 0th spec affected region
# 128: (Added) Unknown artifacts
# 256: Filled CCDs gaps by gfextract/apall
##############
# DQ=8 is considered as good pixels after correction, which is individually identified for Sci.(the above section) and Flat (this section).
# DQ=32 is considered as good pixels. Only correct for amp 4/5/8/9 of Flat. 
# DQ=64(4) is corrected for Flat, and just keep the value in Sci. 
# DQ>=128, i.e., DQ=128 or 256 are considered as bad pixels, which should not be used in the final output. 

prefix = 'xrg'
# fit DQ=32 (edge) and mask DQ=64 (saturated/0th-spec) for Flat
# Do not add DQ for edges at this step to avoid bad fit.
for wave in CentWaves :
    print(wave)
    for obs in Flat[wave] : 
        print(prefix+obs)
        iraf.delete( 'p'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'p'+prefix+obs+'.fits' )
        for amp in range(1, 12+1) : 
            iraf.imreplace( 'p'+prefix+obs+'.fits[VAR,'+str(amp)+']', 0 )
        # linear fit the out edges:
        for (amps, col0, col1) in zip([[1,5,9], [4,8,12]], ['1:4', '253:256'], ['5:16', '241:252']) : 
            for amp in amps : 
                iraf.imdelete( 'tmp_edge_sci' )        
                iraf.fit1d( prefix+obs+'.fits[SCI,'+str(amp)+']', 'tmp_edge_sci', "fit", axis=1, sample=col1, naverage=1, func='chebyshev', order=2, low=3, high=2.3, niterate=5, grow=0, interactive='no' )
                iraf.imcopy( 'tmp_edge_sci.fits['+col0+',*]', 'p'+prefix+obs+'.fits[SCI,'+str(amp)+',overwrite]['+col0+',*]' )
                iraf.imcalc( prefix+obs+'.fits[SCI,'+str(amp)+'], p'+prefix+obs+'.fits[SCI,'+str(amp)+']', 'p'+prefix+obs+'.fits[VAR,'+str(amp)+',overwrite]', 'im1-im2', verbose='no' )
        # join amplifers to identify and fit-fix / growdq the 0th-spec affected region: 
        iraf.delete( 'p'+prefix+obs+'_join.fits' )
        tmp = joinamp( 'p'+prefix+obs, ['SCI', 'VAR', 'DQ'], num=3 ) # output = _join
        # extending saturated (DQ=4) region by 4 pixels (radius=4), change 4 to 64
        iraf.imreplace( 'p'+prefix+obs+'_join.fits[DQ,9]', 64, upper=4, lower=4, radius=4)

prefix = 'pxrg'
# mask DQ=8 and fit-correct outliers for Flat
for wave in CentWaves :
    print(wave)
    for obs in Flat[wave] : 
        print(prefix+obs)
#         if obs!='N20210516S0055' : continue
        # read gain and read noise from header
        tmp = fits.open( prefix+obs+'.fits' )
        RDNOISE = np.array( [ tmp[amp].header['RDNOISE'] for amp in range(1, 12+1) ] )
        GAIN = np.array( [ tmp[amp].header['GAIN'] for amp in range(1, 12+1) ] )
        # calculate mean (used below) and std, and check the ratio of std/mean
        print('GAIN mean, std, std/mean', np.mean(GAIN), np.std(GAIN), np.std(GAIN)/np.mean(GAIN))
        print('RDNOISE mean, std, std/mean', np.mean(RDNOISE), np.std(RDNOISE), np.std(RDNOISE)/np.mean(RDNOISE))
        for amp in [1,5,9] : 
            iraf.imdelete( 'tmp_fit1d' )
            iraf.fit1d( prefix+obs+'_join.fits[SCI,'+str(amp)+']', 'tmp_fit1d', "fit", bpm=prefix+obs+'_join.fits[DQ,'+str(amp)+']', axis=1, sample="*", naverage=1, func='spline3', order=128, low=3, high=2.3, niterate=100, grow=0, interactive='no' )
            # use func='spline3', order=128 instead of func='chebyshev', order=40
            # amp1/5/9 in _join keep amp1/5/9's headers, which do not have LTV1=-32
            # (even amp have, e.g., amp2), so do not need to wcscopy(DQ, SCI)
            ########################################
            iraf.imdelete( 'tmp_residual' )
            iraf.imarith( prefix+obs+'_join.fits[SCI,'+str(amp)+']', '-', 'tmp_fit1d', 'tmp_residual' )
            iraf.imdelete( 'tmp_resmed5' )
            iraf.median( 'tmp_residual', 'tmp_resmed5', 5, 5, zlor=INDEF, zhir=INDEF, verbose='no')
            iraf.imarith( 'tmp_resmed5', "+", 'tmp_fit1d', 'tmp_resmed5')
            ########################################
            iraf.imdelete( 'tmp_noise' )
            iraf.imcalc('tmp_resmed5', 'tmp_noise', "sqrt(im1*"+str(np.mean(GAIN))+"+"+str(np.mean(RDNOISE))+"**2)/"+str(np.mean(GAIN)), verbose='no')
            iraf.imdelete( 'tmp_sigma' )
            iraf.imarith( 'tmp_residual', '/', 'tmp_noise', 'tmp_sigma' )
            ########################################
            sigclip=10
            print( 'positive replacement in flat (res/noise > '+str(sigclip)+' sigma): ' )
            iraf.imstat( 'tmp_sigma', fields="npix", lower=sigclip, format='no')
            print( 'negative replacement in flat: (res/noise < -'+str(sigclip)+' sigma)' )
            iraf.imstat( 'tmp_sigma', fields="npix", upper=-sigclip, format='no')
            iraf.imdelete( 'tmp_flat_bpm' )
            iraf.imcopy( 'tmp_sigma', 'tmp_flat_bpm', verbose='no')
            iraf.imreplace( 'tmp_flat_bpm', 0, lower=-sigclip, upper=sigclip )
            iraf.imreplace( 'tmp_flat_bpm', 1, lower=sigclip )
            iraf.imreplace( 'tmp_flat_bpm', 1, upper=-sigclip )
            ########################################
            iraf.imcalc( prefix+obs+'_join.fits[DQ,'+str(amp)+']'+', tmp_flat_bpm', prefix+obs+'_join.fits[DQ,'+str(amp)+', overwrite]', 'im1+8*im2', verbose='no' )
            # now CR-like regions have DQ=8
            # add DQ=32 for edges of amp 1/5/9 and 4/8/12
            iraf.imdelete( 'tmp_edge_dq' )
            iraf.imcalc( prefix+obs+'_join.fits[DQ,'+str(amp)+'][1:4,*]', 'tmp_edge_dq', 'im1+32', verbose='no' )
            iraf.imcopy( 'tmp_edge_dq.fits[DQ]', prefix+obs+'_join.fits[DQ,'+str(amp)+', overwrite][1:4,*]' )
            iraf.imdelete( 'tmp_edge_dq' )
            iraf.imcalc( prefix+obs+'_join.fits[DQ,'+str(amp)+'][1021:1024,*]', 'tmp_edge_dq', 'im1+32', verbose='no' )
            iraf.imcopy( 'tmp_edge_dq.fits[DQ]', prefix+obs+'_join.fits[DQ,'+str(amp)+', overwrite][1021:1024,*]' )
            ########################################
            iraf.imcalc( prefix+obs+'_join.fits[VAR,'+str(amp)+'],'+prefix+obs+'_join.fits[DQ,'+str(amp)+'], tmp_residual', prefix+obs+'_join.fits[VAR,'+str(amp)+', overwrite]', 'im1+(im2>0)*im3', verbose='no' )
            iraf.imcalc( prefix+obs+'_join.fits[SCI,'+str(amp)+'],'+prefix+obs+'_join.fits[DQ,'+str(amp)+'], tmp_residual', prefix+obs+'_join.fits[SCI,'+str(amp)+', overwrite]', 'im1-(im2>0)*im3', verbose='no' )

# copy the corrected block back to original MEF format file
for wave in CentWaves :
    print(wave)
    for obs in Flat[wave] : 
        print(prefix+obs)
        for amp in [1,5,9] : 
            for i in range(4) : 
                xpos = str(1+256*i)+':'+str(256+256*i)
                iraf.imcopy( prefix+obs+'_join.fits[SCI,'+str(amp)+']['+xpos+',*]', prefix+obs+'.fits[SCI,'+str(amp+i)+'][*,*]' )
                iraf.imcopy( prefix+obs+'_join.fits[VAR,'+str(amp)+']['+xpos+',*]', prefix+obs+'.fits[VAR,'+str(amp+i)+'][*,*]', verbose='no' )
                iraf.imcopy( prefix+obs+'_join.fits[DQ,'+str(amp)+']['+xpos+',*]', prefix+obs+'.fits[DQ,'+str(amp+i)+'][*,*]', verbose='no' )

# delete tmp files
os.system( ' rm -v tmp* ' )

#############
# check results:
# just re-run if showing errors
prefix = 'pxrg'
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] :
        print(prefix+obs)
        iraf.delete( prefix+obs+'_join.fits' )
        tmp = joinamp( prefix+obs, ['SCI', 'VAR', 'DQ'] )
        iraf.display( prefix+obs+'_join.fits[SCI]', 1 )
        iraf.display( prefix+obs+'_join.fits[VAR]', 2, zscale='no', zrange='no', z1=-5, z2=1e4)
        raw_input( ) # pause
        iraf.display( prefix+obs+'_join.fits[DQ]', 2, zscale='no', zrange='no', z1=-5, z2=100)
        raw_input( ) # pause
#############

##############################

prefix = 'xrg'
# mask DQ for Sci. 
for wave in CentWaves :
    print(wave)
    for obs in Sci[wave] : 
        print(prefix+obs)
        iraf.delete( 'p'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'p'+prefix+obs+'.fits')
        # mask DQ=32 (edge) for Sci
        for (amps, col) in zip([[1, 5, 9], [4, 8, 12]], ['1:4', '253:256']) : 
            for amp in amps : 
                iraf.imdelete( 'tmp_edge_dq' )
                iraf.imcalc( 'p'+prefix+obs+'.fits[DQ,'+str(amp)+']['+col+',*]', 'tmp_edge_dq', 'im1+32', verbose='no' )
                iraf.imcopy( 'tmp_edge_dq.fits[DQ]', 'p'+prefix+obs+'.fits[DQ,'+str(amp)+',overwrite]['+col+',*]' )
        # copy DQ=64 (0th-spec) for Sci. using the DQ of corresponding Flat
        for amp in range(1, 12+1) : 
            iraf.imcalc( 'p'+prefix+nearest_file(obs, Flat[wave])+'.fits[DQ,'+str(amp)+']'+', '+'p'+prefix+obs+'.fits[DQ,'+str(amp)+']', 'p'+prefix+obs+'.fits[DQ,'+str(amp)+', overwrite]', '(im1 >= 64)*64+im2', verbose='no')

# check if an additional artifact line exists in long-wavelength direction of 0th order.
# Artifact line show features even extending to the top/bottom (non-exposure) pixels, 
# which should be marked in DQ (review using ds9, frame_panel: tile). 
# If not marked (i.e., additional artifact line), then mark manually with implot. 
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] :
        print(prefix+obs)
        iraf.delete( 'p'+prefix+obs+'_join.fits' )
        tmp = joinamp( 'p'+prefix+obs, ['SCI', 'DQ'] )
        iraf.display( 'p'+prefix+obs+'_join.fits[SCI]', 1 )
        iraf.display( 'p'+prefix+obs+'_join.fits[DQ]', 2, zscale='no', zrange='no', z1=-5, z2=70)
        raw_input( ) # pause # ztrans='log'

for wave in CentWaves:
    print(wave)
    for obs in Sci[wave]: 
        print(obs)
        for col in ['10','11']: 
            iraf.implot( 'p'+prefix+obs+'.fits[SCI,'+col+']', 77 )
            iraf.implot( 'p'+prefix+obs+'.fits[DQ,'+col+']', 77 )
# l for line; c for column; a-a for average

# # Additional artifact line DQ=1(weak) or 128(severe)
# if DQ=128, which should be ejected when combining
# if DQ=1, can be used
###################################################
#### Depends on different data ####
# n1
iraf.imreplace( 'p'+prefix+'N20210306S0051.fits[DQ, 10][256,*]', 1)
iraf.imreplace( 'p'+prefix+'N20210306S0051.fits[DQ, 11][1:7,*]', 1)
iraf.imreplace( 'p'+prefix+'N20210306S0052.fits[DQ, 10][256,*]', 1)
iraf.imreplace( 'p'+prefix+'N20210306S0052.fits[DQ, 11][1:7,*]', 1)
iraf.imreplace( 'p'+prefix+'N20210306S0053.fits[DQ, 10][256,*]', 1)
iraf.imreplace( 'p'+prefix+'N20210306S0053.fits[DQ, 11][1:7,*]', 1)
###################################################

#############
# check results:
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] :
        print(prefix+obs)
        iraf.delete( 'p'+prefix+obs+'_join.fits' )
        tmp = joinamp( 'p'+prefix+obs, ['SCI', 'VAR', 'DQ'] )
        iraf.display( 'p'+prefix+obs+'_join.fits[SCI]', 1 )
        iraf.display( 'p'+prefix+obs+'_join.fits[VAR]', 2, zscale='no', zrange='no', z1=-5, z2=100)
        raw_input( ) # pause
        iraf.display( 'p'+prefix+obs+'_join.fits[DQ]', 2, zscale='no', zrange='no', z1=-5, z2=70)
        raw_input( ) # pause # ztrans='log'

##############################################################
##############################################################

###### Remove scattered light for Flat and Sci. ######

prefix = 'pxrg'

# join amplifers
# just re-run if showing error
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] + Sci[wave] :
        print(prefix+obs)
        iraf.delete( prefix+obs+'_join.fits' )
        tmp = joinamp( prefix+obs, ['SCI', 'VAR', 'DQ'], num=1 ) # output = _join
        
for wave in CentWaves:
    print(wave)
#     if wave=='620' : continue
#     for obs in Flat[wave] :
    for obs in Flat[wave] + Sci[wave] :
        print(prefix+obs)
        # copy DQ
        iraf.delete( 'tmp_sca_bpm.fits')
        iraf.imcopy( prefix+obs+'_join.fits[DQ]', 'tmp_sca_bpm' )
        # change DQ so that gaps have DQ=0
        iraf.imreplace( 'tmp_sca_bpm', 1024)
        blkmask = np.genfromtxt( 'blkmask_rg'+nearest_file(obs, Flat[wave]), skip_header = 0, names = [ 'col0', 'col1', 'line0', 'line1' ], dtype = 'i8, i8, i8, i8' ) 
        for ind in range(0, len(blkmask)):
            iraf.imreplace( 'tmp_sca_bpm.fits[*, '+str(blkmask['line0'][ind])+':'+str(blkmask['line1'][ind])+']', 0)
        iraf.delete( 'tmp_sca_spline.fits')
        # SCI header have LVT2=-48, copy it to bpm otherwise a shift exists in fit1d.
        iraf.wcscopy( 'tmp_sca_bpm', prefix+obs+'_join.fits[SCI]')
        # fit scattered light with spline function
        iraf.fit1d( prefix+obs+'_join.fits[SCI]', 'tmp_sca_spline', 'fit', bpm='tmp_sca_bpm', axis=2, sample="*", naverage=1, func='spline3', order=6, low=3, high=2.3, niterate=5, grow=0, interactive='no' ) 
        # for spline3 a minimum number of data group is 4, so the order is 16/(4-1)=5.3~6
        # 16 is the number of gaps, -1 to consider the adjecant spline pieces
        # subtract scattered light and record it to VAR
        iraf.delete( 'b'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'b'+prefix+obs+'.fits' )
        for amp in range(1,12+1) : 
        	xpos = str(1+256*(amp-1))+':'+str(256+256*(amp-1))
        	iraf.imarith( prefix+obs+'.fits[SCI,'+str(amp)+']', '-', 'tmp_sca_spline.fits['+xpos+',*]', 'b'+prefix+obs+'.fits[SCI,'+str(amp)+',overwrite]' )
        	iraf.imarith( prefix+obs+'.fits[VAR,'+str(amp)+']', '+', 'tmp_sca_spline.fits['+xpos+',*]', 'b'+prefix+obs+'.fits[VAR,'+str(amp)+',overwrite]' )

# delete tmp files
os.system( ' rm -v tmp* ' )

#############
# check results:
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] :
        print(prefix+obs)
        iraf.delete( 'b'+prefix+obs+'_join.fits' )
        tmp = joinamp( 'b'+prefix+obs, ['SCI', 'VAR', 'DQ'] )
        iraf.display( 'b'+prefix+obs+'_join.fits[SCI]', 1 )
        iraf.display( 'b'+prefix+obs+'_join.fits[VAR]', 2, zscale='no', zrange='no', z1=5e2, z2=50e2)
        raw_input( ) # pause # z2=100 for Sci and 1e4 for Flat
    for obs in Sci[wave] :
        print(prefix+obs)
        iraf.delete( 'b'+prefix+obs+'_join.fits' )
        tmp = joinamp( 'b'+prefix+obs, ['SCI', 'VAR', 'DQ'] )
        iraf.display( 'b'+prefix+obs+'_join.fits[SCI]', 1 )
        iraf.display( 'b'+prefix+obs+'_join.fits[VAR]', 2, zscale='no', zrange='no', z1=1, z2=20)
        raw_input( ) # pause # z2=100 for Sci and 1e4 for Flat

##############################################################
##############################################################

###### QE correction for Flat and Std/Sci ######

prefix = 'bpxrg'
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] + Sci[wave] : 
        print(prefix+obs)
        refArc = nearest_file(obs,Arc[wave])
        refFlat = nearest_file(obs,Flat[wave])
        iraf.hedit(  'erg'+refArc+'.fits[0]', 'MASKNAME', 'MDF_'+refFlat, add=yes, verify=no) 
        # gqecorr require the same MDF filename. Note that the MDFs need to be the same actually for each wave. 
        iraf.imdelete( 'q'+prefix+obs )
        iraf.imdelete( 'qecorr'+'erg'+refArc )
        iraf.gqecorr( prefix+obs, \
            refimage='erg'+refArc, \
            fl_correct='yes', fl_vardq='yes', verbose='yes', logfile = 'gqecorr.log')
        iraf.hedit(  'erg'+refArc+'.fits[0]', 'MASKNAME', 'MDF_'+nearest_file(refArc,Flat[wave]), add=yes, verify=no) 
        # recover the original MDF filename

##############################################################
##############################################################

###### Extract Flat and Std/Sci ######

prefix = 'qbpxrg'
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] + Sci[wave] : 
        print(prefix+obs)
        os.system( 'rm -v database/ap'+'e'+prefix+obs+'*' )
        iraf.imdelete( 'e'+prefix+obs )
        iraf.gfextract( prefix+obs, \
            reference='erg'+nearest_file(obs,Flat[wave]), \
            recenter='no', trace='no', weights='none', fl_vardq='yes', fl_fulldq='yes')
            # use fl_fulldq='yes' to retain  DQ information correctly

############
# check results:
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] + Sci[wave] : 
        print(prefix+obs)
        iraf.display('e'+prefix+obs+'.fits[SCI, 1]', 1)
        iraf.display('e'+prefix+obs+'.fits[VAR, 1]', 2, zscale='no', zrange='no', z1=5, z2=100)
        raw_input( ) # pause
        iraf.display('e'+prefix+obs+'.fits[DQ, 1]', 2, zscale='no', zrange='no', z1=-1, z2=64)
        raw_input( ) # pause
        iraf.display('e'+prefix+obs+'.fits[SCI, 2]', 1)
        iraf.display('e'+prefix+obs+'.fits[VAR, 2]', 2, zscale='no', zrange='no', z1=5, z2=100)
        raw_input( ) # pause
        iraf.display('e'+prefix+obs+'.fits[DQ, 2]', 2, zscale='no', zrange='no', z1=-1, z2=64)
        raw_input( ) # pause

# check the cross-matching between fiber ID and its locaiton
for wave in CentWaves :
    print(wave)
    for obs in Sci[wave] :
        print(obs)
        iraf.gfdisplay( 'e'+prefix+obs, 1, z1=15e4, z2=60e4, ztrans='linear' )

##############################################################
##############################################################

###### Fix pixels in CCD gaps for Sci. (skip Flat) ######

# It is suggested that gfextract uses fixpix (linear) to fill those bad pixels, 
# which works well for Flat but sometimes bad for Sci. if the adjacent regions
# are CR effected or with bright sky lines. 
# Therefore, here re-fix these bad pixels with fit1d in axis=1, a low highrej with 
# a high niter are used to force the fitting line near to the continum component. 
# Temparaly do not fix saturated regions (DQ=64).

prefix='eqbpxrg'
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] : 
        print(prefix+obs)
        iraf.imdelete( 'p'+prefix+obs )
        iraf.gemfix( prefix+obs, 'p'+prefix+obs, bitmask=256, axis=1, method="fit1d", order=10, low_reject=3, high_reject=2.3, niterate=100, grow=0, fl_inter='no' ) 
        iraf.tcopy( prefix+obs+'.fits[MDF]', 'p'+prefix+obs+'.fits[MDF, overwrite]', verbose='yes')
        # recover MDF changed by gemfix bug
    for obs in Flat[wave]:
        print(prefix+obs)
        iraf.delete( 'p'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'p'+prefix+obs+'.fits' )

############
# check results: 
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] : 
        print(obs) 
        iraf.display( prefix+obs+'.fits[SCI, 1]', 1 )
        iraf.display( 'p'+prefix+obs+'.fits[SCI, 1]', 2 )
        raw_input( ) # pause
        iraf.display( prefix+obs+'.fits[SCI, 2]', 1 )
        iraf.display( 'p'+prefix+obs+'.fits[SCI, 2]', 2 )
        raw_input( ) # pause

##############################################################
##############################################################

###### Rectify the spectra for Flat and Std/Sci ######

prefix='peqbpxrg'
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] + Sci[wave] : 
        print(prefix+obs)
        iraf.imdelete( 't'+prefix+obs)
        # iraf.gftransform( prefix+obs, wavtraname='erg'+nearest_file(obs,Arc[wave]), fl_vardq='no')
        iraf.gftransform( prefix+obs, wavtraname='erg'+nearest_file(obs,Arc[wave]), w1=4525, w2=10804, dw=3.9, fl_vardq='yes')

###### Wavelength Alignment using Sky lines ######

prefix='tpeqbpxrg'
expr_sky = "XINST > 10"
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] : 
        print(prefix+obs)
        iraf.imdelete( 'a'+prefix+obs+'_sky' )
        iraf.gfapsum( prefix+obs, outimages='a'+prefix+obs+'_sky' , expr=expr_sky, combine='sum', fl_inter='no')

# interactively:
iraf.delete( 'rvLog.txt' )
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] : 
        print(prefix+obs)
        iraf.delete( 'database/id'+'a'+prefix+obs+'_sky' )
        rv.observatory='gemini-north'
        rv.rvidlines ( 'a'+prefix+obs+'_sky[SCI]',  \
            coordlist=obs_dir+'auxiliary/wavecalib/skylines_ex.txt', nsum=1, maxfeatures=10, \
            ftype='emission', fwidth=10, cradius=10, threshold=7, minsep=5, \
            logfile='rvLog.txt', autowrite=yes)
# l to identify the night sky lines; the peaks in 5500--9000AA range
# f to fit for the mean velocity
# q to write the parameters to the log file and quit the task
# the output z is average of (λ_obs-λ_true)/λ_true, 
 
 # check the residuals of calibration:
iraf.type( 'rvLog.txt' )
###################################################
#### Depends on different data ####
# the (λ_obs-λ_true) is within [-2,2] with rms of 1--2, which is smaller than dw=3.9
# so just keep the above rectified image without further correction.
# the rms of velocity is 50--80 km/s, can be considered as velocity accuracy.
###################################################

############
# check results: 
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] : 
        print(prefix+obs)
        iraf.display( prefix+obs+'.fits[SCI]', 1)
        iraf.display( prefix+obs+'.fits[VAR]', 2, zscale='no', zrange='no', z1=8, z2=80)
        raw_input( ) # pause
        iraf.display( prefix+obs+'.fits[DQ]', 2, zscale='no', zrange='no', z1=-1, z2=72)
        raw_input( ) # pause
    for obs in Flat[wave] : 
        print(prefix+obs)
        iraf.display( prefix+obs+'.fits[SCI]', 1)
        iraf.display( prefix+obs+'.fits[VAR]', 2, zscale='no', zrange='no', z1=800, z2=16000)
        raw_input( ) # pause
        iraf.display( prefix+obs+'.fits[DQ]', 2, zscale='no', zrange='no', z1=-1, z2=72)
        raw_input( ) # pause

##############################################################
##############################################################

###### Correct for 0th spectra for Flat ######

prefix='tpeqbpxrg'

# since the red-slit has gaps (only one CentWave), here use reduced standards flat
# as reference and skip the normalizing of red-slit for the same grism/filter set. 
# reference = '../auxiliary/tpeqbpxrgS20BFT210_flat'
reference_1 = obs_dir+'auxiliary/flat/ptpeqbpxrgFeige34_flat'
reference_2 = obs_dir+'auxiliary/flat/ptpeqbpxrgKF08T3_flat'
os.system( 'rm -v tmp_ref*' )
iraf.imcombine( reference_1+'.fits[SCI], '+reference_2+'.fits[SCI]', 'tmp_ref.fits', \
    sigma='tmp_ref_sigma.fits', \
    reject='none', lthreshold=1e-6, blank=0 ) 

for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] : 
        print(prefix+obs)
        # calculate ratio of Flat on refered one
        iraf.delete( prefix+obs+'_calib.fits' )
        iraf.imarith( prefix+obs+'.fits[SCI]', '/', 'tmp_ref.fits', prefix+obs+'_calib.fits' )
        ########
        # fit 0th spec region using adjecent regions: 
        iraf.delete( 'tmp_calib_fit1d.fits' )
        iraf.fit1d( prefix+obs+'_calib.fits', 'tmp_calib_fit1d.fits', "fit", axis=1, sample='430:440, 540:550', naverage=1, func='chebyshev', order=2, low=3, high=2.3, niterate=10, grow=0, interactive='no' )
        iraf.imcopy( 'tmp_calib_fit1d.fits[440:540, 749:1494]', prefix+obs+'_calib.fits[440:540, 749:1494]' )
        ########
        # fill gaps: 
        iraf.delete( 'tmp_bpm.fits' )
        max_allowed_dq = 127
        iraf.imcalc( prefix+obs+'[SCI], '+prefix+obs+'[DQ]', 'tmp_bpm.fits', '1-(im1>1e-6)*(im2<='+str(max_allowed_dq)+')' )
        #iraf.imreplace( 'tmp_bpm.fits', 1, lower=1, radius=2 )
        iraf.fixpix( prefix+obs+'_calib.fits', 'tmp_bpm.fits', linterp=1 )
        ########
        # correct the original Flat file
        iraf.delete( 'p'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'p'+prefix+obs+'.fits' )
        iraf.imarith( prefix+obs+'_calib.fits', '*', 'tmp_ref.fits', 'p'+prefix+obs+'.fits[SCI, 1,overwrite]' )
        iraf.delete( 'tmp_var.fits' )
        iraf.imarith( 'p'+prefix+obs+'.fits[SCI]', '-', prefix+obs+'.fits[SCI]', 'tmp_var.fits' )
        iraf.imarith( 'p'+prefix+obs+'.fits[VAR]', '+', 'tmp_var.fits', 'p'+prefix+obs+'.fits[VAR,overwrite]' )

##################
# check corrected Flat:
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] : 
        print(prefix+obs)
        iraf.display( 'p'+prefix+obs+'.fits[SCI]', 1 )
        # iraf.display( 'p'+prefix+obs+'.fits[VAR]', 2 )
        iraf.display( prefix+obs+'.fits[SCI]', 2 )
        raw_input( ) # pause

##############################################################
##############################################################

###### Response function from 0th-corrected Flat ######

# In order to keep consistency with Feige34, i.e., spectral calibration, 
# here use the sum spectrum from the red slit region of combined Feige34 Flat. 

prefix='ptpeqbpxrg'
resp_spec=obs_dir+'auxiliary/flat/aptpeqbpxrgFeige34_flat_red'

iraf.delete( 'tmp_mean_1d.fits' )
iraf.imcopy( resp_spec+'.fits[SCI]', 'tmp_mean_1d.fits' )
iraf.delete( 'tmp_mean_2d.fits' )
iraf.imstack( 'tmp_mean_1d.fits', 'tmp_mean_2d.fits' )
iraf.blkrep( 'tmp_mean_2d.fits', 'tmp_mean_2d.fits', 1, fits.open( prefix+Flat[CentWaves[0]][0]+'.fits' )[2].data.shape[0] )

# derive response images
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] : 
        print(prefix+obs)
        iraf.delete( prefix+obs+'_resp.fits' )
        iraf.imarith( prefix+obs+'.fits[SCI]', '/', 'tmp_mean_2d.fits', prefix+obs+'_resp.fits[SCI]' )

# check the image:
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] : 
        print(prefix+obs)
        iraf.display( prefix+obs+'[SCI]', 1)
        iraf.display( prefix+obs+'_resp[SCI]', 2, zscale='no', zrange='no', z1=0.75, z2=1.75)
        raw_input( ) # pause

##############################################################
##############################################################

###### Fiber response correction for Flat (only for check) and Sci ######

prefix='tpeqbpxrg'

for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] + Sci[wave]  : 
        print(prefix+obs)
        iraf.delete( 'i'+prefix+obs+'.fits')
        iraf.copy( prefix+obs+'.fits', 'i'+prefix+obs+'.fits')
        resp_image = 'p'+prefix+nearest_file(obs,Flat[wave])+'_resp.fits'
        iraf.imarith( 'i'+prefix+obs+'.fits[SCI]', '/', resp_image+'[SCI]', 'i'+prefix+obs+'.fits[SCI,overwrite]' )

# check the image:
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] + Sci[wave]  : 
        print(prefix+obs)
        iraf.display( 'i'+prefix+obs+'.fits[SCI]', 1)#, zscale='no', zrange='no', z1=0.75, z2=1.75)
        raw_input( ) # pause

##############################################################
##############################################################

###### Clean remaining CR (mainly in y-axis) ######

prefix='itpeqbpxrg'
# os.system( 'rm -v lacos*' )

line0s = [1, 749] # start position of red/blue slit
line1s = [748, 1494] # last position of red/blue slit
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] :
        print(obs)
        iraf.delete( 'x'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'x'+prefix+obs+'.fits' )
        ########################################
        for (line0, line1) in zip(line0s, line1s): 
            print(line0, line1)
            ypos = str(line0)+':'+str(line1)
            # use modified lacos_spec to obtain clean SCI and mask
            iraf.imdelete( 'tmpclean1' )
            iraf.imdelete( 'tmpmask1' )
            iraf.lacos_spec_mod( prefix+obs+'.fits[SCI][*,'+ypos+']', 'tmpclean1', 'tmpmask1', gain=1, readn=0, xorder=20, yorder=10, sigclip=10, sigfrac=0.5, objlim=0.5, niter=4, verbose='no', keeptmp='no')
            # use sigclip=10 to protect bright sky and galaxy (e.g, Hα) lines (default sigclip=4.5)
            # use objlim=0.5 to select CR near to bright sky and galaxy lines (default objlim=2)
            ########################################
            # caiculate the CR flux and recored into VAR
            iraf.imdelete( 'tmpvar1' )
            iraf.imarith( prefix+obs+'.fits[SCI][*,'+ypos+']', '-', 'tmpclean1', 'tmpvar1' )
            # scale tmpmask1 to DQ=8 and add into original DQ 
            iraf.imreplace( 'tmpmask1', 8, lower=1 )
            iraf.imarith( prefix+obs+'.fits[DQ][*,'+ypos+']', '+', 'tmpmask1', 'tmpmask1' )
            ########################################
            # copy the corrected block back to original MEF format file
            iraf.imcopy( 'tmpclean1.fits[*,*]', 'x'+prefix+obs+'.fits[SCI][*,'+ypos+']' )
            iraf.imcopy( 'tmpvar1.fits[*,*]', 'x'+prefix+obs+'.fits[VAR][*,'+ypos+']', verbose='no' )
            iraf.imcopy( 'tmpmask1.fits[*,*]', 'x'+prefix+obs+'.fits[DQ][*,'+ypos+']', verbose='no' )

# just copy file for Flat
prefix='itpeqbpxrg'
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave]:
        print(obs)
        iraf.delete( 'x'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'x'+prefix+obs+'.fits' )

#########################
# NOTE: be careful to check if strong galaxy lines (e.g., Hα, [OIII]) are identified as CR. 
# If so, use a larger objlim, e.g., 3--5 
#########################

#############
# check results:
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] :
        print(obs)
        iraf.display( prefix+obs+'.fits[SCI]', 1 )
        iraf.display( 'x'+prefix+obs+'.fits[VAR]', 2, zscale='no', zrange='no', z1=-5, z2=100)
        raw_input( ) # pause
        iraf.display( 'x'+prefix+obs+'.fits[SCI]', 1 )
        raw_input( ) # pause
        
##############################################################
##############################################################

###### Spectrophotometric calibration for Flat (only for check) and Sci ######

prefix='xitpeqbpxrg'
sens_image = obs_dir+'auxiliary/sensitivity/asitpeqbpxrgFeige34_sens_v3'
extinct_file = obs_dir+'auxiliary/sensitivity/Ext_cont_line.dat'

for wave in CentWaves:
    print(wave)
    for (filelist, fl_ext) in zip( [Flat[wave], Sci[wave]], ['no', 'yes'] ):
        for obs in filelist:
            print(prefix+obs)
            iraf.imdelete( 'c'+prefix+obs)
            iraf.gscalibrate( prefix+obs, \
                sfunction=sens_image, \
                obs='Gemini-North', \
                extinction=extinct_file, fl_ext=fl_ext, fl_vardq='yes')
                # use extinction curve including telluric absorption

#############
# check results:
prefix='xitpeqbpxrg'
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] + Sci[wave]  : 
        print(prefix+obs)
        iraf.display( 'c'+prefix+obs+'.fits[SCI]', 1)#, zscale='no', zrange='no', z1=0.75, z2=1.75)
        iraf.display( 'c'+prefix+obs+'.fits[DQ]', 2, zscale='no', zrange='no', z1=0, z2=72)
        raw_input( ) # pause

##############################################################
##############################################################

###### Creation of model of unknown artifact features around 8000A ######

# only run for J0941 to obtain a template
# skip this section if use a reference unk_model

##############################################################
##############################################################

###### Correction of unknown artifact features around 8000A ######

prefix='cxitpeqbpxrg'

unk_Dir_ref = obs_dir+'auxiliary/artifactmod/'
unk_mod_ref = 'cxitpeqbpxrgJ0941_n1_comb_unk_mod.fits'
iraf.delete( 'tmp_unk_mod.fits' )
iraf.imcopy( unk_Dir_ref+unk_mod_ref, 'tmp_unk_mod.fits' )

# determine APIDs_unk_nogal (artifact region without galaxy emission)
# and APIDs_ref (region with similar background level of artifact)
# create tmp cube 
obs = Sci[CentWaves[0]][0]      
iraf.imdelete( 'tmp4gfcube')
iraf.copy( prefix+obs+'.fits', 'tmp4gfcube.fits' )
# add some non-zero values into the starting and ending wavelegnth positions, 
# otherwise gfcube will not convert the full wavelegnth range.
iraf.imreplace( 'tmp4gfcube.fits[SCI][1, *]', -1e6 )
iraf.imreplace( 'tmp4gfcube.fits[SCI][1611, *]', -1e6 )
iraf.imdelete( 'tmp_cube')
iraf.gfcube( 'tmp4gfcube', outimage='tmp_cube', \
    ssample=0.05, bitmask=0, \
    fl_atmdisp='yes', fl_flux='yes', fl_var='no', fl_dq='no') # logfile='gfcubeLog.txt'
# set the values at the starting and ending wavelegnth positions back to zero. 
iraf.imreplace( 'tmp_cube.fits[SCI][*, *, 1]', 0 )
iraf.imreplace( 'tmp_cube.fits[SCI][*, *, 1611]', 0 )
########
# determine inds_spec_unk_all and stack range (whether in CCD gaps)
iraf.implot(prefix+obs+'.fits[SCI]', 1072) # BLOCK 6-25
# create tmp artifact stack image
iraf.delete( 'tmp_wavestack_unk.fits' )
iraf.imcombine( 'tmp_cube.fits[SCI][*,*,800:930]', 'tmp_wavestack_unk.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' ) 
iraf.display( 'tmp_wavestack_unk.fits', 2, zscale='no', zrange='no', z1=0.015, z2=0.025) # 0.02, 0.03 for n2
# determine BLOCKs with above snapshots comparing with hex grids
# check the APIDs with the BLOCKs
tread(prefix+obs+'.fits[MDF]') 
# view MDF, ^D + goto column; ^D + quit to quit
###################################################
#### Depends on different data ####
APIDs_unk_nogal = range(1093, 1096+1) + range(1147, 1149+1) 
# BLOCK: 6-4 -- 6-1, 5-50 -- 5-47 (5-49 bad)
APIDs_ref = range(1190, 1193+1) + range(1196, 1199+1) # BLOCK: 5-5 -- 5-2, 4-49 -- 4-46
############
inds_spec_unk_all = range(800, 930+1)
###################################################

for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] : 
        print(obs)
        # obtain an average spectrum in refered APIDs
        iraf.delete( 'tmp_unk_ref1d.fits' )
        iraf.imcombine( prefix+obs+'.fits[SCI][*,'+str(min(APIDs_ref))+':'+str(max(APIDs_ref))+']', 'tmp_unk_ref1d.fits', reject='none', project='yes' ) 
        iraf.delete( 'tmp_unk_ref2d.fits' )
        iraf.imstack( 'tmp_unk_ref1d.fits', 'tmp_unk_ref2d.fits' )
        iraf.blkrep( 'tmp_unk_ref2d.fits', 'tmp_unk_ref2d.fits', 1, 1494 )
        # reduce the averaged reference spectrum from data
        iraf.delete( 'tmp_unk_residual.fits' )
        iraf.imarith( prefix+obs+'.fits[SCI]', '-', 'tmp_unk_ref2d.fits', 'tmp_unk_residual.fits' )
        # copy bpm
        iraf.delete( 'tmp_unk_bpm.fits' )
        iraf.imcopy( prefix+obs+'.fits[DQ]', 'tmp_unk_bpm.fits' )
        ##########
        iraf.delete( 'tmp_unk_res_gap9.fits' )
        iraf.imcalc( 'tmp_unk_residual.fits, tmp_unk_bpm.fits', 'tmp_unk_res_gap9.fits', 'im1*(im2<32)-9999*(im2>=32)' )
        iraf.delete( 'tmp_unk_res_ws.fits' )
        iraf.imtranspose( 'tmp_unk_res_gap9.fits', 'tmp_unk_res_gap9.fits' )
        iraf.imcombine( 'tmp_unk_res_gap9.fits[*,'+str(min(inds_spec_unk_all))+':'+str(max(inds_spec_unk_all))+']', 'tmp_unk_res_ws.fits', reject='none', lthreshold=-9e3, blank=0, project='yes' )
        ##########
        iraf.delete( 'tmp_unk_mod_gap9.fits' )
        iraf.imcalc( 'tmp_unk_mod.fits, tmp_unk_bpm.fits', 'tmp_unk_mod_gap9.fits', 'im1*(im2<32)-9999*(im2>=32)' )
        iraf.delete( 'tmp_unk_mod_ws.fits' )
        iraf.imtranspose( 'tmp_unk_mod_gap9.fits', 'tmp_unk_mod_gap9.fits' )
        iraf.imcombine( 'tmp_unk_mod_gap9.fits[*,'+str(min(inds_spec_unk_all))+':'+str(max(inds_spec_unk_all))+']', 'tmp_unk_mod_ws.fits', reject='none', lthreshold=-9e3, blank=0, project='yes' )
        ##########
        # obtain the ratio of non-galaxy fibers over unk_model
        iraf.delete( 'tmp_unk_ratio.fits' )
        iraf.imarith( 'tmp_unk_res_ws.fits', '/', 'tmp_unk_mod_ws.fits', 'tmp_unk_ratio.fits' )
        iraf.delete( 'tmp_unk_nogal.fits' )
        iraf.imarith( 'tmp_unk_ratio.fits', '*', 0, 'tmp_unk_nogal.fits' )
        for APID in APIDs_unk_nogal : 
            iraf.imreplace( 'tmp_unk_nogal.fits[0]['+str(APID)+']', 1 )
            # add extvar[0] otherwise with header error
        iraf.imarith( 'tmp_unk_ratio.fits', '*', 'tmp_unk_nogal.fits', 'tmp_unk_ratio.fits' )
        tmp_data = fits.open( 'tmp_unk_ratio.fits' )[0].data
        unk_ratio = np.mean(tmp_data[ tmp_data > 0 ])
        print(unk_ratio)
        # created the modeled unk_spectra for subtraction
        iraf.delete( prefix+obs+'_unk_mod.fits' )
        iraf.imcalc( 'tmp_unk_mod.fits', prefix+obs+'_unk_mod.fits', 'im1*'+str(unk_ratio)+'' )
        # subtract the modeled unk_spectra
        iraf.delete( 'p'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'p'+prefix+obs+'.fits' )
        iraf.imarith( prefix+obs+'.fits[SCI]', '-', prefix+obs+'_unk_mod.fits' , 'p'+prefix+obs+'.fits[SCI, overwrite]' )
        iraf.imarith( prefix+obs+'.fits[VAR]', '+', prefix+obs+'_unk_mod.fits' , 'p'+prefix+obs+'.fits[VAR, overwrite]' )

###################################################
#### Depends on different data ####
# mannually correct
obs = 'N20210306S0052'
unk_ratio = 0.8 # from mean of 0051/0.82 and 0053/0.78
# created the modeled unk_spectra for subtraction
iraf.delete( prefix+obs+'_unk_mod.fits' )
iraf.imcalc( 'tmp_unk_mod.fits', prefix+obs+'_unk_mod.fits', 'im1*'+str(unk_ratio)+'' )
# subtract the modeled unk_spectra
iraf.delete( 'p'+prefix+obs+'.fits' )
iraf.copy( prefix+obs+'.fits', 'p'+prefix+obs+'.fits' )
iraf.imarith( prefix+obs+'.fits[SCI]', '-', prefix+obs+'_unk_mod.fits' , 'p'+prefix+obs+'.fits[SCI, overwrite]' )
iraf.imarith( prefix+obs+'.fits[VAR]', '+', prefix+obs+'_unk_mod.fits' , 'p'+prefix+obs+'.fits[VAR, overwrite]' )
###################################################

# just copy Flat (not affected, possibly due to short exposure)
for wave in CentWaves:
    print(wave)
    for obs in Flat[wave] : 
        print(obs)
        iraf.delete( 'p'+prefix+obs+'.fits' )
        iraf.copy( prefix+obs+'.fits', 'p'+prefix+obs+'.fits' )

#############
# check results:
prefix='cxitpeqbpxrg'
obs = Sci[CentWaves[0]][1]      
iraf.imdelete( 'tmp4gfcube')
iraf.copy( 'p'+prefix+obs+'.fits', 'tmp4gfcube.fits' )
iraf.imreplace( 'tmp4gfcube.fits[SCI][1, *]', -1e6 )
iraf.imreplace( 'tmp4gfcube.fits[SCI][1611, *]', -1e6 )
iraf.imdelete( 'tmp_cube')
iraf.gfcube( 'tmp4gfcube', outimage='tmp_cube', \
    ssample=0.05, bitmask=0, \
    fl_atmdisp='yes', fl_flux='yes', fl_var='no', fl_dq='no') 
iraf.imreplace( 'tmp_cube.fits[SCI][*, *, 1]', 0 )
iraf.imreplace( 'tmp_cube.fits[SCI][*, *, 1611]', 0 )
iraf.delete( 'tmp_wavestack_unk.fits' )
iraf.imcombine( 'tmp_cube.fits[SCI][*,*,800:930]', 'tmp_wavestack_unk.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' ) 
iraf.display( 'tmp_wavestack_unk.fits', 1, zscale='no', zrange='no', z1=0.015, z2=0.025)
# for wave in CentWaves:
#     print(wave)
#     for obs in Sci[wave]  : 
#         print(prefix+obs)
#         iraf.display( 'p'+prefix+obs+'.fits[SCI]', 1)
#         iraf.display( 'p'+prefix+obs+'.fits[VAR]', 2)
#         raw_input( ) # pause
#         iraf.display( 'p'+prefix+obs+'.fits[DQ]', 2, zscale='no', zrange='no', z1=0, z2=72)
#         raw_input( ) # pause

##############################################################
##############################################################

###### Preparation for combining Sci. Cubes ######

prefix='pcxitpeqbpxrg'

# set the max allowed DQ
# exclude DQ=128, 256
# keep DQ=8, 32, 64(keep 0th region)
max_allowed_dq = 127 

for wave in CentWaves:
    print(wave)
    for obs in Sci[wave] : 
        if (obs != 'N20210306S0052') : continue
        print(obs)
        iraf.imdelete( prefix+obs+'_norm' )
        iraf.copy( prefix+obs+'.fits', prefix+obs+'_norm.fits' )
        iraf.imreplace( prefix+obs+'_norm.fits[DQ]', 1024, lower=1+max_allowed_dq )
        iraf.imreplace( prefix+obs+'_norm.fits[DQ]', 1, lower=0, upper=max_allowed_dq )
        iraf.imreplace( prefix+obs+'_norm.fits[DQ]', 0, lower=1024, upper=1024 )
        # now good:bad pixels in _norm,fits have DQ=1:0
        # use corresponds Flat to identify non-empty regions
        iraf.imcopy( prefix+nearest_file(obs,Flat[wave])+'.fits[SCI][*,*]', prefix+obs+'_norm.fits[SCI][*,*]'  )
        iraf.imreplace( prefix+obs+'_norm.fits[SCI]', 1, lower=1e-6 )
        iraf.imarith( prefix+obs+'_norm.fits[DQ]', '*', prefix+obs+'_norm.fits[SCI]', prefix+obs+'_norm.fits[DQ, overwrite]', pixtype='ushort' )
        iraf.imreplace( prefix+obs+'_norm.fits[DQ]', 1, lower=1 )
        iraf.imreplace( prefix+obs+'_norm.fits[DQ]', 0, upper=1-1e-6 )
        # now empty-region pixels in _norm have DQ=0
        # rewrite to SCI for combine
        iraf.imcopy( prefix+obs+'_norm.fits[DQ]', prefix+obs+'_norm.fits[SCI][*,*]' )
        iraf.imdelete( prefix+obs+'_gap0' )
        iraf.copy( prefix+obs+'.fits', prefix+obs+'_gap0.fits' )
        iraf.imarith( prefix+obs+'_gap0.fits[SCI]', '*', prefix+obs+'_norm.fits[DQ]', prefix+obs+'_gap0.fits[SCI, overwrite]', pixtype='real' )

# if with error: ERROR (502, "floating point invalid operation"), just repeat 
#############
# check results:
prefix='pcxitpeqbpxrg'
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave]  : 
        print(obs)
        iraf.display( prefix+obs+'_gap0.fits[SCI]', 1)#, zscale='no', zrange='no', z1=0.75, z2=1.75)
        iraf.display( prefix+obs+'_gap0.fits[DQ]', 2, zscale='no', zrange='no', z1=2, z2=36)
        raw_input( ) # pause
        iraf.display( prefix+obs+'_norm.fits[SCI]', 2, zscale='no', zrange='no', z1=-0.05, z2=1.05)
        raw_input( ) # pause

##############################################################
##############################################################

###### Create the cubes of Sci. ######

prefix='pcxitpeqbpxrg'
for wave in CentWaves :
    print(wave)
    for obs in Sci[wave] : 
        for suffix in ['_gap0', '_norm'] : 
            if (obs != 'N20210306S0052') : continue
            print(prefix+obs+suffix)
            iraf.imdelete( 'tmp4gfcube')
            iraf.copy( prefix+obs+suffix+'.fits', 'tmp4gfcube.fits' )
            # add some non-zero values into the starting and ending wavelegnth positions, 
            # otherwise gfcube will not convert the full wavelegnth range.
            iraf.imreplace( 'tmp4gfcube.fits[SCI][1, *]', -1e6 )
            iraf.imreplace( 'tmp4gfcube.fits[SCI][1611, *]', -1e6 )
            iraf.imdelete( 'd'+prefix+obs+suffix)
            iraf.gfcube( 'tmp4gfcube', outimage='d'+prefix+obs+suffix, \
                ssample=0.05, bitmask=0, \
                fl_atmdisp='yes', fl_flux='yes', fl_var='no', fl_dq='no', \
                logfile='gfcubeLog.txt')
            # set the values at the starting and ending wavelegnth positions back to zero. 
            iraf.imreplace( 'd'+prefix+obs+suffix+'.fits[SCI][*, *, 1]', 0 )
            iraf.imreplace( 'd'+prefix+obs+suffix+'.fits[SCI][*, *, 1611]', 0 )

# normalize _norm to unit of 1 and exclude the regions with _norm < 0.9
iraf.minmax( 'd'+prefix+Sci[CentWaves[0]][0]+'_norm.fits[SCI][*,*, 200]')
iraf.minmax( 'd'+prefix+Sci[CentWaves[0]][0]+'_norm.fits[SCI][*,*, 1200]')
nmax = 28.86750030517578 # depends on ssample of gfcube
nlim = 0.9 # set a strict threshold to remove edge pixels with poor quality
for wave in CentWaves :
    print(wave)
    for obs in Sci[wave] : 
        if (obs != 'N20210306S0052') : continue
        print(prefix+obs)
        # normalize _norm to unit of 1
        iraf.imarith( 'd'+prefix+obs+'_norm.fits[SCI]', '/', nmax, 'd'+prefix+obs+'_norm.fits[SCI,overwrite]' )
        # exclude _gap0 with _norm < nlim and re-normalize with _norm (oscillattion at edges)
        iraf.imcalc( 'd'+prefix+obs+'_gap0.fits[SCI],'+'d'+prefix+obs+'_norm.fits[SCI]', \
            'd'+prefix+obs+'_gap0.fits[SCI,overwrite]', \
            'im1*(im2>='+str(nlim)+')/im2+0*(im2<'+str(nlim)+')', verbose='no' )
        # re-normalize _norm to bi-value 0 and 1
        iraf.imcalc( 'd'+prefix+obs+'_norm.fits[SCI]', 'd'+prefix+obs+'_norm.fits[SCI,overwrite]', \
            '1*(im1>='+str(nlim)+')+0*(im1<'+str(nlim)+')', verbose='no' )

# set valid pixels fraction and remove invalid points:
flim_y = 0.75 # for fitting of sky subtraction
flim_x = 0.25 # exclude to few valid points; should < 0.50 to allow single slit coverage
for wave in CentWaves : 
    print(wave)
    for obs in Sci[wave] : 
        if (obs != 'N20210306S0052') : continue
        print(prefix+obs)
        norm = fits.open( 'd'+prefix+obs+'_norm.fits' )[1].data
        for ind_spec in range(1, norm.shape[0]+1) : #using iraf list index
            #if ind_spec!=350: continue
            # print(ind_spec)
            nsum_xy = np.sum( norm[ind_spec-1, :, :] ) # note python list index
            isinvalid_spec = nsum_xy < (flim_y*norm.shape[1]*flim_x*norm.shape[2])
            # print( nsum_xy/norm.shape[1]/norm.shape[2], 'invalid' if isinvalid_spec else 'normal' )
            if ( isinvalid_spec ) : 
                iraf.imreplace( 'd'+prefix+obs+'_norm.fits[SCI][*, *, '+str(ind_spec)+']', 0 )
                continue
            ######
            nsum_x = np.sum( norm[ind_spec-1, :, :], 1 )
            ind_invalid_y = 1+np.where( nsum_x < (flim_x*norm.shape[2]) )[0]
            if ( ind_invalid_y.shape[0] > 0 ) : 
                ind_invalid_yp, ind_invalid_yn = ind_invalid_y.copy(), ind_invalid_y.copy()
                ind_invalid_yp[0:-1], ind_invalid_yp[-1] = ind_invalid_y[1:], 1e6
                ind_invalid_yn[0], ind_invalid_yn[1:] = -1e6, ind_invalid_y[0:-1]
                ind_invalid_ystart = ind_invalid_y[ np.where( (ind_invalid_y-ind_invalid_yn) != 1 )[0] ]
                ind_invalid_ylast = ind_invalid_y[ np.where( (ind_invalid_y-ind_invalid_yp) != -1 )[0] ]
                for ind_ygroup in range( ind_invalid_ystart.shape[0] ) : 
                    ind_y = str(ind_invalid_ystart[ind_ygroup])+':'+str(ind_invalid_ylast[ind_ygroup])
                    # print(ind_y)
                    iraf.imreplace( 'd'+prefix+obs+'_norm.fits[SCI][*, '+ind_y+', '+str(ind_spec)+']', 0 )
            ######
            nsum_y = np.sum( norm[ind_spec-1, :, : ], 0 )
            ind_invalid_x = 1+np.where( nsum_y < (flim_y*norm.shape[1]) )[0]            
            if ( ind_invalid_x.shape[0] > 0 ) : 
                ind_invalid_xp, ind_invalid_xn = ind_invalid_x.copy(), ind_invalid_x.copy()
                ind_invalid_xp[0:-1], ind_invalid_xp[-1] = ind_invalid_x[1:], 1e6
                ind_invalid_xn[0], ind_invalid_xn[1:] = -1e6, ind_invalid_x[0:-1]
                ind_invalid_xstart = ind_invalid_x[ np.where( (ind_invalid_x-ind_invalid_xn) != 1 )[0] ]
                ind_invalid_xlast = ind_invalid_x[ np.where( (ind_invalid_x-ind_invalid_xp) != -1 )[0] ]
                for ind_xgroup in range( ind_invalid_xstart.shape[0] ) : 
                    ind_x = str(ind_invalid_xstart[ind_xgroup])+':'+str(ind_invalid_xlast[ind_xgroup])
                    # print(ind_x)
                    iraf.imreplace( 'd'+prefix+obs+'_norm.fits[SCI]['+ind_x+', *, '+str(ind_spec)+']', 0 )
        ############
        # remove invalid points in _gap0
        iraf.imcalc( 'd'+prefix+obs+'_gap0.fits[SCI],'+'d'+prefix+obs+'_norm.fits[SCI]', \
        'd'+prefix+obs+'_gap0.fits[SCI,overwrite]', 'im1*im2', verbose='no' )      

#############
# check results:
prefix='pcxitpeqbpxrg'
for wave in CentWaves:
    print(wave)
    for obs in Sci[wave]  : 
        print(obs)
        iraf.display( 'd'+prefix+obs+'_gap0.fits[SCI][*,*, 865]', 1)#, zscale='no', zrange='no', z1=0, z2=0.03) #375
        iraf.display( 'd'+prefix+obs+'_norm.fits[SCI][*,*, 865]', 2, zscale='no', zrange='no', z1=-0.5, z2=30)
        raw_input( ) # pause

# prefix='pcxitpeqbpxrg'
# obs = 'N20210306S0053'
# for ind_spec in range(820, 900)  : 
#     print(ind_spec)
#     iraf.display( 'd'+prefix+obs+'_gap0.fits[SCI][*,*, '+str(ind_spec)+']', 1)#, zscale='no', zrange='no', z1=0, z2=0.03) #375
#     raw_input( ) # pause
# 
# iraf.display( 'd'+prefix+obs+'_gap0.fits[SCI][28,*,*]', 1, zscale='no', zrange='no', z1=0, z2=0.03) #375

##############################################################
##############################################################

###### Create wavelength stacking images ######

prefix='dpcxitpeqbpxrg'
obses = [ obs for sublist in [Sci[wave] for wave in CentWaves] for obs in sublist ]
xsize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[2]
ysize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[1]
wsize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[0]

nlim = 0.75
for obs in obses : 
    if (obs != 'N20210306S0052') : continue
    iraf.delete( prefix+obs+'_framecover.fits' )
    iraf.blkavg( prefix+obs+'_norm.fits[SCI]', prefix+obs+'_framecover.fits', xsize, ysize, 1, option='average' )
    iraf.blkrep( prefix+obs+'_framecover.fits', prefix+obs+'_framecover.fits', xsize, ysize, 1 )
    iraf.delete( prefix+obs+'_gap9.fits' )
    iraf.imcalc( prefix+obs+'_gap0.fits[SCI], '+prefix+obs+'_norm.fits[SCI], '+prefix+obs+'_framecover.fits', prefix+obs+'_gap9.fits', 'im1*(im2>='+str(nlim)+')*(im3>='+str(nlim)+')-9999*((im2<'+str(nlim)+')+(im3<'+str(nlim)+'))', verbose='no' )
    # avoid ind_spec=480:520 region with 0th-order contaminations
    # avoid ind_spec=720:1000 region with unknown contaminations
    iraf.delete( prefix+obs+'_wavestack_p1.fits' )
    iraf.imcombine( prefix+obs+'_gap9.fits[SCI][*,*,250:450]', prefix+obs+'_wavestack_p1.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' ) #lsigma=3, hsigma=3, nkeep=1, 
    iraf.delete( prefix+obs+'_wavestack_p2.fits' )
    iraf.imcombine( prefix+obs+'_gap9.fits[SCI][*,*,520:720]', prefix+obs+'_wavestack_p2.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' )
    iraf.delete( prefix+obs+'_wavestack_p3.fits' )
    iraf.imcombine( prefix+obs+'_gap9.fits[SCI][*,*,1000:1200]', prefix+obs+'_wavestack_p3.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' )
    iraf.delete( prefix+obs+'_wavestack.fits' )
    iraf.imcombine( prefix+obs+'_wavestack_p1.fits, '+prefix+obs+'_wavestack_p2.fits, '+prefix+obs+'_wavestack_p3.fits', prefix+obs+'_wavestack.fits', reject='none', lthreshold=-9e3, blank=-9999 )
    # stack ind_spec=815:915 region with unknown contaminations (check the correction)
    iraf.delete( prefix+obs+'_wavestack_unk.fits' )
    iraf.imcombine( prefix+obs+'_gap9.fits[SCI][*,*,820:900]', prefix+obs+'_wavestack_unk.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' ) #lsigma=3, hsigma=3, nkeep=1, 

#############
# check results:
# stack ind_spec=815:915 region with unknown contaminations (check the correction)
prefix='dpcxitpeqbpxrg'
for obs in obses :  
    iraf.display( prefix+obs+'_wavestack.fits', 1)#, zscale='no', zrange='no', z1=0.008, z2=0.02)
    iraf.display( prefix+obs+'_wavestack_p1.fits', 2)
    iraf.display( prefix+obs+'_wavestack_p2.fits', 3)
    iraf.display( prefix+obs+'_wavestack_p3.fits', 4)
    raw_input( ) # pause
    iraf.display( prefix+obs+'_wavestack_unk.fits', 2)
    raw_input( ) # pause
    
##############################################################
##############################################################

###### Subtract sky lines ######

prefix='dpcxitpeqbpxrg'
obses = [ obs for sublist in [Sci[wave] for wave in CentWaves] for obs in sublist ]
xsize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[2]
ysize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[1]
wsize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[0]
    
# use sky line free image to set fitting areas:
for obs in obses : 
    print(prefix+obs)
    iraf.display( prefix+obs+'_wavestack.fits', 1, zscale='no', zrange='no', z1=0.015, z2=0.03)
    raw_input( ) # pause
# select the rough flux level of bound of galaxy region
###################################################
#### Depends on different data ####
# thresholds = [0.0175, 0.0185, 0.0195] # n1
thresholds = [0.0225, 0.0220, 0.0220] # n2
###################################################
# check tmp_bpm using the above threshold
for (obs, threshold) in zip(obses, thresholds) : 
    iraf.delete( prefix+obs+'_skymod_bpm_init.fits' )
    iraf.imcalc( prefix+obs+'_wavestack.fits', prefix+obs+'_skymod_bpm_init.fits' , '(im1>'+str(threshold)+')' )
    iraf.display( prefix+obs+'_wavestack.fits', 1, zscale='no', zrange='no', z1=0.015, z2=0.03)
    iraf.display( prefix+obs+'_skymod_bpm_init.fits', 3)
    raw_input( ) # pause
##############
# create a unified bpm via overlapping
# NOTE that this assumes that all images have the same PSF and no large alignment shifts.
iraf.delete( 'tmp_bpm_overlap.fits' )
iraf.imcombine( ','.join( [prefix+obs+'_skymod_bpm_init.fits' for obs in obses ] ), 'tmp_bpm_overlap.fits', reject='none' )
iraf.display( 'tmp_bpm_overlap.fits', 3 )
iraf.imreplace( 'tmp_bpm_overlap.fits', 0, upper=0.75 )
iraf.imreplace( 'tmp_bpm_overlap.fits', 1, lower=0.75, radius=6 ) # 2 for n2
# modify grow-radius to allow tolarence due to variance and shifts
iraf.display( 'tmp_bpm_overlap.fits', 4 )
###################################################
#### Depends on different data ####
iraf.imreplace( 'tmp_bpm_overlap.fits[*,1:7]', 0 ) # n2
###################################################
##############
# copy bpm to cube format adding _norm constraint
for (obs, threshold) in zip(obses, thresholds) : 
    iraf.delete( 'tmp_bpm_cube.fits' )
    iraf.imstack( 'tmp_bpm_overlap.fits', 'tmp_bpm_cube.fits' )
    iraf.blkrep( 'tmp_bpm_cube.fits', 'tmp_bpm_cube.fits', 1, 1, wsize )
    iraf.delete( prefix+obs+'_skymod_bpm.fits' )
    iraf.imcalc( prefix+obs+'_norm.fits[SCI], '+'tmp_bpm_cube.fits', prefix+obs+'_skymod_bpm.fits[SCI]', '1-(1-im2)*(im1>0.75)' )
    # put the input with EXTNAME at the begining, 
    # otherwise the output could not get EXTNAME.
##############
# # check bpm map:
# for obs in obses : 
#     iraf.display( prefix+obs+'_skymod_bpm.fits[SCI][*,*,120]', 1)
#     iraf.display( prefix+obs+'_skymod_bpm.fits[SCI][*,*,1200]', 2)
#     raw_input( ) # pause

############################
############################

# if REPEAT, start from here

# NOTE that fit1d sometimes fails with a non-rectangular bpm map, i.e., in a given wavelength range,
# therefore hereafter perform fit1d in y-pos per given x-pos and wavelength. 
# (Note that the bpm map should also be rectangular in y-pos vs. wavelength plane.)
# The fit is along y-axis, since sky line mainly shows features along y-axis. 
# Fit in y-axis is good enough, do not consider a 2nd fit in x-axis or smoothing sky map. 

fit_order = 2 # linear
for obs in obses : 
    iraf.delete( prefix+obs+'_skymod_y.fits' )
    iraf.imarith( prefix+obs+'_gap0.fits[SCI]', '*', 0, prefix+obs+'_skymod_y.fits[SCI]' )
    for ind_x in range(1, xsize+1) : 
        #ind_x = 30
        print(ind_x)
        iraf.delete( 'tmp_prefit_yw.fits' )
        iraf.imcopy( prefix+obs+'_gap0.fits[SCI]['+str(ind_x)+',*,*]', 'tmp_prefit_yw.fits', verbose='no' )
        iraf.hedit( 'tmp_prefit_yw.fits', 'WCSDIM, WAXMAP01', 0, delete=yes, verify=no, show=no) 
        iraf.delete( 'tmp_bpm_yw.fits' )
        iraf.imcopy( prefix+obs+'_skymod_bpm.fits[SCI]['+str(ind_x)+',*,*]', 'tmp_bpm_yw.fits', verbose='no' )
        iraf.hedit( 'tmp_bpm_yw.fits', 'WCSDIM, WAXMAP01', 0, delete=yes, verify=no, show=no)
        # delete header keywords to avoid fit1d bug: MWCS: dimension mismatch (mw_gltermd) 
        ######################################
        list_bpm_0 = fits.open( 'tmp_bpm_yw.fits' )[0].data
        # estimate valid (bpm=0) pixels in each (x,spec) slice at both ends 
        nstart = int(ysize/2) - np.sum( list_bpm_0[:,0:int(ysize/2)], 1 )
        nlast = int(ysize/2) - np.sum( list_bpm_0[:,(ysize-int(ysize/2)):ysize], 1 )
        # for given x-pos, calculate the starting and ending wavelength with the same valid bpm
        list_bpm_p, list_bpm_n = list_bpm_0.copy(), list_bpm_0.copy()
        list_bpm_p[0:-1,:], list_bpm_p[-1,:] = list_bpm_0[1:,:], 1e6
        list_bpm_n[0,:], list_bpm_n[1:,:] = -1e6, list_bpm_0[0:-1,:]
        ind_spec_start = 1+np.where( (np.sum(np.abs(list_bpm_0-list_bpm_n),1) != 0) * (nstart >= 3) * (nlast >= 3) )[0]
        ind_spec_last = 1+np.where( (np.sum(np.abs(list_bpm_0-list_bpm_p),1) != 0) * (nstart >= 3) * (nlast >= 3) )[0]
        #############################################
        if (ind_spec_start.shape[0] == 0) : 
            print( 'no enough sky bkg' )
            continue
        for ind_spec_group in range( ind_spec_start.shape[0] ) : 
            ind_spec = str(ind_spec_start[ind_spec_group])+':'+str(ind_spec_last[ind_spec_group])
            print(ind_spec)
            iraf.delete( 'tmp_prefit_y.fits' )
            iraf.imcopy( 'tmp_prefit_yw.fits[*,'+ind_spec+']', 'tmp_prefit_y.fits', verbose='no' )
            iraf.delete( 'tmp_bpm_y.fits' )
            iraf.imcopy( 'tmp_bpm_yw.fits[*,'+ind_spec+']', 'tmp_bpm_y.fits', verbose='no' )
            iraf.delete( 'tmp_skymod_y.fits' )
            iraf.fit1d( 'tmp_prefit_y.fits', 'tmp_skymod_y.fits', bpm='tmp_bpm_y.fits', \
                type='fit', interactive='no', axis=1, naverage=1, \
                function='chebyshev', order=fit_order, low_reject=3, high_reject=2.3, niterate=10, grow=0 )
            #iraf.imcalc( 'tmp_skymod_y.fits,'+prefix+obs+'_norm.fits[SCI]['+str(ind_x)+',*,'+ind_spec+']', 'tmp_skymod_y.fits', 'im1*(im2>0.75)', verbose='no' )
            iraf.imcopy( 'tmp_skymod_y.fits', prefix+obs+'_skymod_y.fits[SCI]['+str(ind_x)+',*,'+ind_spec+']', verbose='no' )

# subtract sky map from _gap0
for obs in obses : 
    iraf.delete( 's'+prefix+obs+'_gap0.fits' )
    #iraf.copy( prefix+obs+'_gap0.fits', 's'+prefix+obs+'_gap0.fits' )
    iraf.imarith( prefix+obs+'_gap0.fits[SCI]', '-', prefix+obs+'_skymod_y.fits[SCI]', 's'+prefix+obs+'_gap0.fits[SCI]' )
    # copy _norm
    iraf.delete( 's'+prefix+obs+'_norm.fits' )
    iraf.copy( prefix+obs+'_norm.fits', 's'+prefix+obs+'_norm.fits' )

# stack along wavelength 
nlim = 0.75
for obs in obses : 
    iraf.delete( 's'+prefix+obs+'_gap9.fits' )
    iraf.imcalc( 's'+prefix+obs+'_gap0.fits[SCI], '+prefix+obs+'_norm.fits[SCI], '+prefix+obs+'_framecover.fits', 's'+prefix+obs+'_gap9.fits', 'im1*(im2>='+str(nlim)+')*(im3>='+str(nlim)+')-9999*((im2<'+str(nlim)+')+(im3<'+str(nlim)+'))', verbose='no' )
    # avoid ind_spec=450:550 region with 0th-order contaminations
    # avoid ind_spec=760:960 region with unknown contaminations
    iraf.delete( 's'+prefix+obs+'_wavestack_p1.fits' )
    iraf.imcombine( 's'+prefix+obs+'_gap9.fits[SCI][*,*,250:450]', 's'+prefix+obs+'_wavestack_p1.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' ) #lsigma=3, hsigma=3, nkeep=1, 
    iraf.delete( 's'+prefix+obs+'_wavestack_p2.fits' )
    iraf.imcombine( 's'+prefix+obs+'_gap9.fits[SCI][*,*,525:725]', 's'+prefix+obs+'_wavestack_p2.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' )
    iraf.delete( 's'+prefix+obs+'_wavestack_p3.fits' )
    iraf.imcombine( 's'+prefix+obs+'_gap9.fits[SCI][*,*,1000:1200]', 's'+prefix+obs+'_wavestack_p3.fits', reject='none', lthreshold=-9e3, blank=-9999, project='yes' )
    iraf.delete( 's'+prefix+obs+'_wavestack.fits' )
    iraf.imcombine(  's'+prefix+obs+'_wavestack_p1.fits, '+'s'+prefix+obs+'_wavestack_p2.fits, '+'s'+prefix+obs+'_wavestack_p3.fits', 's'+prefix+obs+'_wavestack.fits', reject='none', lthreshold=-9e3, blank=-9999 )
    iraf.display( 's'+prefix+obs+'_wavestack.fits', 1)
    ############
    iraf.delete( 'tmp_bpm.fits' )
    iraf.imcombine( prefix+obs+'_skymod_bpm.fits[SCI][*,*,250:1200]', 'tmp_bpm.fits', reject='none', project='yes' )
    iraf.display( 'tmp_bpm.fits', 3)
    iraf.delete( 'tmp_imstat.fits' )
    iraf.imcalc( 's'+prefix+obs+'_wavestack.fits, tmp_bpm.fits', 'tmp_imstat.fits', 'im1*(im2<=0.1)-9999*(im2>0.1)', verbose='no' )
    iraf.display( 'tmp_imstat.fits', 2, zscale='no', zrange='no', z1=-0.001, z2=0.01)
    tmp = fits.open( 'tmp_imstat.fits' )
    immean = np.mean(tmp[0].data[ tmp[0].data > -1 ])
    imstddev = np.std(tmp[0].data[ tmp[0].data > -1 ])
    threshold = immean+imstddev*3
    print( 'immean, imstddev, threshold:', immean, imstddev, threshold )
    # re-generate bpm using the above value and 3-sigma threshold
    iraf.delete( prefix+obs+'_skymod_bpm_init.fits' )
    iraf.imcalc( 's'+prefix+obs+'_wavestack.fits', prefix+obs+'_skymod_bpm_init.fits', '(im1>'+str(immean+imstddev*3)+')', verbose='no' )
    iraf.display( prefix+obs+'_skymod_bpm_init.fits', 4)
    # raw_input( ) # pause

# create a unified bpm via overlapping
# NOTE that this assumes that all images have the same PSF.
iraf.display( 's'+prefix+obs+'_wavestack.fits', 1)
iraf.display( 'tmp_imstat.fits', 2, zscale='no', zrange='no', z1=-0.001, z2=0.01)
iraf.delete( 'tmp_bpm_overlap.fits' )
iraf.imcombine( ','.join( [prefix+obs+'_skymod_bpm_init.fits' for obs in obses ] ), 'tmp_bpm_overlap.fits', reject='none' )
iraf.display( 'tmp_bpm_overlap.fits', 3 )
iraf.imreplace( 'tmp_bpm_overlap.fits', 0, upper=0.75 )
iraf.imreplace( 'tmp_bpm_overlap.fits', 1, lower=0.75, radius=4 ) # 3 for n2
# modify grow-radius to allow tolarence due to variance and shifts
# # mannually correct edges
# iraf.imreplace( 'tmp_bpm_overlap.fits[1,*]', 0 )
# iraf.imreplace( 'tmp_bpm_overlap.fits['+str(xsize)+',*]', 0 )
# iraf.imreplace( 'tmp_bpm_overlap.fits[*, 1]', 0 )
# iraf.imreplace( 'tmp_bpm_overlap.fits[*, '+str(ysize)+']', 0 )
iraf.display( 'tmp_bpm_overlap.fits', 4 )

# if REPEAT above steps, copy bpm to cube format with _norm constraint
for (obs, threshold) in zip(obses, thresholds) : 
    iraf.delete( 'tmp_bpm_cube.fits' )
    iraf.imstack( 'tmp_bpm_overlap.fits', 'tmp_bpm_cube.fits' )
    iraf.blkrep( 'tmp_bpm_cube.fits',  'tmp_bpm_cube.fits', 1, 1, wsize )
    iraf.delete( prefix+obs+'_skymod_bpm.fits' )
    iraf.imcalc( prefix+obs+'_norm.fits[SCI], '+'tmp_bpm_cube.fits', prefix+obs+'_skymod_bpm.fits[SCI]', '1-(1-im2)*(im1>0.75)' )
    # put the input with EXTNAME at the begining, 
    # otherwise the output could not get EXTNAME.

# ##############
# # check bpm map:
# for obs in obses : 
#     iraf.display( prefix+obs+'_skymod_bpm.fits[SCI][*,*,1200]', 4)
#     raw_input( ) # pause

###################################
# repeat sky subtraction using the new bpm for one iteration
###################################

#############
# check results:
for obs in obses : 
    # for ind_spec in range(250, 1450+1) : 
    # for ind_spec in [375,797,950,1105,1232,1346,1353,1391,1404] : 
    for ind_spec in [375,797,860,950,1105,1232,1346,1391,1404] : 
        print(ind_spec)
        tmp = fits.open( prefix+obs+'_gap0.fits' )
        immax = np.max(tmp[1].data[ind_spec,:,:])
        if immax == 0: continue
        iraf.display( prefix+obs+'_gap0.fits[SCI][*,*,'+str(ind_spec)+']', 1, zscale='no', zrange='no', z1=-0.01, z2=immax*0.75)
        iraf.display( prefix+obs+'_skymod_y.fits[SCI][*,*,'+str(ind_spec)+']', 2, zscale='no', zrange='no', z1=-0.01, z2=immax*0.75)
        iraf.display( 's'+prefix+obs+'_gap0.fits[SCI][*,*,'+str(ind_spec)+']', 3, zscale='no', zrange='no', z1=-0.01, z2=immax*0.75)
        iraf.display( prefix+obs+'_skymod_bpm.fits[SCI][*,*,'+str(ind_spec)+']', 4, zscale='no', zrange='no', z1=0, z2=1)
        raw_input( ) # pause

# obs = obses[0]
# ind_x = 65
# iraf.display( prefix+obs+'_gap0.fits[SCI]['+str(ind_x)+',*,*]', 1, zscale='no', zrange='no', z1=-0.01, z2=0.06)
# iraf.display( prefix+obs+'_skymod_y.fits[SCI]['+str(ind_x)+',*,*]', 2, zscale='no', zrange='no', z1=-0.01, z2=0.06)
# iraf.display( 's'+prefix+obs+'_gap0.fits[SCI]['+str(ind_x)+',*,*]', 3, zscale='no', zrange='no', z1=-0.001, z2=0.01)
# iraf.display( prefix+obs+'_skymod_bpm.fits[SCI]['+str(ind_x)+',*,*]', 4)

# spectrum in a central pixel
prefix='dpcxitpeqbpxrg'
for obs in obses : 
    iraf.implot( 's'+prefix+obs+'_gap0.fits[SCI][64, 48, *]')

##############################################################
##############################################################

###### Spatial Alignment ######

prefix='sdpcxitpeqbpxrg'
obses = [ obs for sublist in [Sci[wave] for wave in CentWaves] for obs in sublist ]
xsize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[2]
ysize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[1]
wsize = fits.open( prefix+obses[0]+'_gap0.fits' )[1].data.shape[0]

###################################################
#### Depends on different data ####
# considering n1/n2 observations, use N20210322S0090 as reference
Dir_ref = '../J1126_n2_v72/'
obs_ref = 'N20210322S0090'
###################################################

# search for coordinates of bright object
# iraf.display( prefix+Sci[CentWaves[0]][0]+'_wavestack.fits', 1, zscale='no', zrange='no', z1=0, z2=0.02)
os.system('rm -v imcentroid_coords.dat')
###################################################
#### Depends on different data ####
os.system('echo 57 45 > imcentroid_coords.dat')
###################################################
iraf.type('imcentroid_coords.dat')

# estimate the alignment shifts
iraf.imcentroid( ','.join( [ prefix+obs+'_wavestack' for obs in obses ] ), Dir_ref+prefix+obs_ref+'_wavestack', 'imcentroid_coords.dat', bigbox=25 )
###################################################
#### Depends on different data ####
# n1, 210306
xshifts = [-0.150,  0.355, -0.108]
yshifts = [ 0.730,  0.468, -0.027]
# n2, 210322
xshifts = [ 0.393, -0.455,  0.000]
yshifts = [ 0.044,  0.043,  0.000]
###################################################

# align _gap0 and _norm using the obtained shifts
for (obs, xshift, yshift) in zip(obses, xshifts, yshifts) : 
    print(obs, xshift, yshift)
    for suffix in ['_gap0', '_norm'] : 
        iraf.delete( 'l'+prefix+obs+suffix+'.fits' )
        iraf.copy( prefix+obs+suffix+'.fits', 'l'+prefix+obs+suffix+'.fits' )
        if ( (xshift==0)*(yshift==0) ) : continue
        for ind_spec in range(1, wsize+1): 
            iraf.delete( 'tmp_aligned.fits' )
            iraf.imshift( prefix+obs+suffix+'.fits[SCI][*,*,'+str(ind_spec)+']', 'tmp_aligned.fits', xshift, yshift, interp_type='drizzle[0.1]', boundary_type="constant", constant=0 ) # 
            # use drizzle[0.1] and constant=0 for a better fitting at boundaries.
            # interp_type options: drizzle[0.5], linear, spline3
            # boundary_type options: constant, nearest
            iraf.imcopy( 'tmp_aligned.fits', 'l'+prefix+obs+suffix+'.fits[SCI][*,*,'+str(ind_spec)+']', verbose='no' ) 

# exclude the pixels with valid fraction < 0.75 after alignment
nlim = 0.75
for obs in obses : 
    print(prefix+obs)
    # exclude _gap0 with _norm < nlim and re-normalize with _norm
    iraf.imcalc( 'l'+prefix+obs+'_gap0.fits[SCI],'+'l'+prefix+obs+'_norm.fits[SCI]', \
    'l'+prefix+obs+'_gap0.fits[SCI,overwrite]', \
        'im1*(im2>='+str(nlim)+')/im2+0*(im2<'+str(nlim)+')', verbose='no' )
    # re-normalize _norm to bi-value 0 and 1
    iraf.imcalc( 'l'+prefix+obs+'_norm.fits[SCI]', 'l'+prefix+obs+'_norm.fits[SCI,overwrite]', \
        '1*(im1>='+str(nlim)+')+0*(im1<'+str(nlim)+')', verbose='no' )

#############
# check results:
prefix='sdpcxitpeqbpxrg'
ind_wave, ind_obs, ind_spec = 0, 0, 1050#350
iraf.display( prefix+Sci[CentWaves[ind_wave]][ind_obs]+'_gap0.fits[SCI][*,*,'+str(ind_spec)+']', 1, zscale='no', zrange='no', z1=0, z2=0.01)
iraf.display( 'l'+prefix+Sci[CentWaves[ind_wave]][ind_obs]+'_gap0.fits[SCI][*,*,'+str(ind_spec)+']', 2, zscale='no', zrange='no', z1=0, z2=0.01)
iraf.display( prefix+Sci[CentWaves[ind_wave]][ind_obs]+'_norm.fits[SCI][*,*,'+str(ind_spec)+']', 3, zscale='no', zrange='no', z1=-0.3, z2=1.3)
iraf.display( 'l'+prefix+Sci[CentWaves[ind_wave]][ind_obs]+'_norm.fits[SCI][*,*,'+str(ind_spec)+']', 4, zscale='no', zrange='no', z1=-0.3, z2=1.3)

##############################################################
##############################################################

###### Combine Sci. cubes ###### 

prefix='lsdpcxitpeqbpxrg'
obses = [ obs for sublist in [Sci[wave] for wave in CentWaves] for obs in sublist ]

# change gap (no data) value from 0 to -9999 to use lthreshold in imcombine
nlim = 0.75
for obs in obses : 
    iraf.delete( prefix+obs+'_gap9.fits' )
    iraf.imcalc( prefix+obs+'_gap0.fits[SCI], '+prefix+obs+'_norm.fits[SCI]', prefix+obs+'_gap9.fits', 'im1*(im2>='+str(nlim)+')-9999*(im2<'+str(nlim)+')', verbose='no' )

# combine _gap9
suffix = '_gap9'
os.system( 'rm -v '+prefix+Sci_comb[0]+'*' )
iraf.imcombine( ','.join( [ prefix+obs+suffix+'.fits[SCI]'  for obs in obses ] ), \
    prefix+Sci_comb[0]+'.fits', \
    bpmasks=prefix+Sci_comb[0]+'_bpmasks.fits', \
    rejmask=prefix+Sci_comb[0]+'_rejmask.fits', \
    nrejmasks=prefix+Sci_comb[0]+'_nrejmasks.fits', \
    expmasks=prefix+Sci_comb[0]+'_expmasks.fits', \
    sigma=prefix+Sci_comb[0]+'_sigma.fits', \
    reject='sigclip', lsigma=3, hsigma=3, nkeep=1, lthreshold=-9e3, blank=0 )
# note that sigclip only work for > 3 input pixels
# remove BPM from header otherwise MWCS bug appears when using display
iraf.hedit( prefix+Sci_comb[0]+'.fits', 'BPM', 0, delete=yes, verify=no, show=no)
# note that in the output blank=0, i.e., actual gap value=0
# combine _norm
# for suffix in ['_gap0', '_norm'] : 
suffix = '_norm'
os.system( 'rm -v '+prefix+Sci_comb[0]+suffix+'*' )
iraf.imcombine( ','.join( [ prefix+obs+suffix+'.fits[SCI]'  for obs in obses ] ), \
prefix+Sci_comb[0]+suffix+'.fits', reject='none' )

#############
# check results:
for ind_spec in [375,797,860,950,1105,1174,1232,1346,1353,1391,1404] : 
    print(ind_spec)
    iraf.display( prefix+Sci_comb[0]+'.fits[*,*,'+str(ind_spec)+']', 1, zscale='no', zrange='no', z1=0, z2=0.01)
    iraf.display( prefix+Sci_comb[0]+'_bpmasks.fits[*,*,'+str(ind_spec)+']', 2, zscale='no', zrange='no', z1=-0.1, z2=1.1)
    iraf.display( prefix+Sci_comb[0]+'_nrejmasks.fits[*,*,'+str(ind_spec)+']', 3, zscale='no', zrange='no', z1=-0.1, z2=4.1)
    iraf.display( prefix+Sci_comb[0]+'_sigma.fits[*,*,'+str(ind_spec)+']', 4)
    raw_input( ) # pause

# ind_y = 40
# iraf.display( prefix+Sci_comb[0]+'.fits[*,'+str(ind_y)+',*]', 1, zscale='no', zrange='no', z1=0, z2=0.04)
# iraf.display( prefix+Sci_comb[0]+'_bpmasks.fits[*,'+str(ind_y)+',*]', 2)
# iraf.display( prefix+Sci_comb[0]+'_nrejmasks.fits[*,'+str(ind_y)+',*]', 3)
# iraf.display( prefix+Sci_comb[0]+'_sigma.fits[*,'+str(ind_y)+',*]', 4)

# spectrum in a central pixel
prefix='lsdpcxitpeqbpxrg'
iraf.implot( prefix+Sci_comb[0]+'.fits[57, 45, *]')

###############
# check alignment
prefix='lsdpcxitpeqbpxrg'
obs=Sci_comb[0]
iraf.delete( prefix+obs+'_wavestack.fits' )
iraf.imcombine( prefix+obs+'.fits[*,*,550:1250]', prefix+obs+'_wavestack.fits', reject='none', project='yes' )
# estimate the alignment shifts
iraf.imcentroid( prefix+obs+'_wavestack.fits', Dir_ref+prefix[1:]+obs_ref+'_wavestack', 'imcentroid_coords.dat' )

##############################################################
##############################################################

###### Calculate error spectrum ######

prefix='lsdpcxitpeqbpxrg'
obs=Sci_comb[0]

# create non-source mask
iraf.display( prefix+obs+'_wavestack.fits', 1)
###################################################
#### Depends on different data ####
threshold = 0.0006 # n1
threshold = 0.0006 # n2
###################################################
iraf.imhistogram( prefix+obs+'_wavestack.fits', binwidth=1e-4, z1=-1e-3, z2=threshold*10, logy='no' )
tmp = fits.open( prefix+obs+'_wavestack.fits' )
tmp_data = tmp[0].data
immean = np.median( tmp_data [ tmp_data < threshold ] )
imstddev = np.std( tmp_data [ tmp_data < threshold ] )
threshold = immean+imstddev*2 # use 2-sigma to avoid source contamination
print( 'immean, imstddev, threshold:', immean, imstddev, threshold )
####
iraf.delete( 'tmp_imstat.fits' )
iraf.imcalc( prefix+obs+'_wavestack.fits', 'tmp_imstat.fits', 'im1*(im1<'+str(threshold)+')-9999*(im1>='+str(threshold)+')', verbose='no' )
iraf.display( 'tmp_imstat.fits', 2, zscale='no', zrange='no', z1=-0.001, z2=0.01)
####
iraf.delete( prefix+obs+'_skymask.fits' )
iraf.imcalc( prefix+obs+'_wavestack.fits', prefix+obs+'_skymask.fits', '1*(im1<'+str(threshold)+')+0*(im1>='+str(threshold)+')', verbose='no' )
iraf.imreplace( prefix+obs+'_skymask.fits', 0, upper=0.1, radius=2 )
iraf.display( prefix+obs+'_skymask.fits', 3) #, zscale='no', zrange='no', z1=-0.001, z2=0.01)

# generate error spectrum
iraf.delete( prefix+obs+'_error_mean.fits' )
iraf.imarith( prefix+obs+'.fits', '*', 0, prefix+obs+'_error_mean.fits' )
iraf.delete( prefix+obs+'_error_std.fits' )
iraf.imarith( prefix+obs+'.fits', '*', 0, prefix+obs+'_error_std.fits' )
xmed = np.int(np.round(xsize/2.0))
for ind_spec in range(1, wsize+1) : 
    # print( ind_spec )
    for ind_slit in ['1:'+str(xmed+2), str(xmed+3)+':'+str(xsize)] : # shift edges to blue slit
        iraf.delete( 'tmp_imstat.fits' )
        iraf.imcalc( prefix+obs+'.fits['+ind_slit+',*,'+str(ind_spec)+'], '+prefix+obs+'_bpmasks.fits.pl['+ind_slit+',*,'+str(ind_spec)+'], '+prefix+obs+'_skymask.fits['+ind_slit+', *]', 'tmp_imstat.fits', 'im1*(im2<0.5)*im3-9999*(1-(im2<0.5)*im3)', verbose='no' )
        tmp = fits.open( 'tmp_imstat.fits' )
        mask_valid = tmp[0].data > -1
        if (np.sum(mask_valid) <= 2) : continue
        immean = np.mean(tmp[0].data[ mask_valid ])
        imstddev = np.std(tmp[0].data[ mask_valid ])
        iraf.imreplace( prefix+obs+'_error_mean.fits['+ind_slit+',*,'+str(ind_spec)+']', immean )
        iraf.imreplace( prefix+obs+'_error_std.fits['+ind_slit+',*,'+str(ind_spec)+']', imstddev )

#############
# check results:
iraf.implot( prefix+obs+'_error_mean.fits[70,50,*]' )
iraf.implot( prefix+obs+'_error_mean.fits[71,50,*]' )
iraf.implot( prefix+obs+'_error_std.fits[70,50,*]' )
iraf.implot( prefix+obs+'_error_std.fits[71,50,*]' )

##############################################################
##############################################################

###### Append extend into Cube file ######

prefix='lsdpcxitpeqbpxrg'
obs=Sci_comb[0]

# correct for _error_mean
iraf.delete( 'u'+prefix+obs+'.fits' )
# iraf.imarith( prefix+obs+'.fits', '-', prefix+obs+'_error_mean.fits', 'u'+prefix+obs+'.fits' )
# just copy the spectra before subtracting mean error
iraf.imcopy( prefix+obs+'.fits', 'u'+prefix+obs+'.fits[SCI, 1, append]' )

# add EXT for later analysis
iraf.imcopy( prefix+obs+'_error_std.fits', 'u'+prefix+obs+'.fits[ERRSTD, 2, append]' )
iraf.imcopy( prefix+obs+'_error_mean.fits', 'u'+prefix+obs+'.fits[ERRMN, 3, append]' )
iraf.imcopy( prefix+obs+'_sigma.fits', 'u'+prefix+obs+'.fits[COMBSTD, 4, append]' )
iraf.imcopy( prefix+obs+'_bpmasks.fits.pl', 'u'+prefix+obs+'.fits[BPM, 5, append]' )
iraf.imcopy( prefix+obs+'_nrejmasks.fits.pl', 'u'+prefix+obs+'.fits[NREJ, 6, append]' )
# iraf.hedit( prefix+Sci_comb+'.fits[0]', 'EXTNAME', 'SCI', update=yes, verify=no)

#############
# check results:
for ind_spec in [375,797,950,1105,1232,1346,1353,1391,1404] : 
    print(ind_spec)
    iraf.display( 'u'+prefix+obs+'.fits[SCI][*,*,'+str(ind_spec)+']', 1, zscale='no', zrange='no', z1=0, z2=0.01)
    iraf.display( 'u'+prefix+obs+'.fits[ERRSTD][*,*,'+str(ind_spec)+']', 2, zscale='no', zrange='no', z1=0, z2=0.01)
    iraf.display( 'u'+prefix+obs+'.fits[ERRMN][*,*,'+str(ind_spec)+']', 3, zscale='no', zrange='no', z1=0, z2=0.01)
    iraf.display( 'u'+prefix+obs+'.fits[COMBSTD][*,*,'+str(ind_spec)+']', 4, zscale='no', zrange='no', z1=0, z2=0.01 )
    raw_input( ) # pause

##############################################################
##############################################################

###### edit SCI header to modify WCS to RA/Dec ######

# https://noirlab.edu/science/programs/csdc/usngo/instruments/gmos/faq#ifu Q3

prefix='lsdpcxitpeqbpxrg'
obs=Sci_comb[0]

# search for the center of the galaxy:
iraf.imcentroid( prefix+obs+'_wavestack.fits', prefix+obs+'_wavestack.fits', 'imcentroid_coords.dat' )
###################################################
#### Depends on different data ####
CRPIX1, CRPIX2 = 57.052, 45.406 # n1
CRPIX1, CRPIX2 = 57.474, 45.282 # n2
###################################################

# read parameters from headers
tmp = fits.open( prefix[2:]+Sci[CentWaves[0]][0]+'_gap0.fits' )
RA_notuse, DEC_notuse, PA = tmp[0].header['RA'], tmp[0].header['DEC'], tmp[0].header['PA']
NAXIS1, NAXIS2, resample = tmp[1].header['NAXIS1'], tmp[1].header['NAXIS2'], abs(tmp[1].header['CD1_1'])/3600 
# resample = arcsec per pixel
CD3_3 = tmp[1].header['CD3_3']

# use alma observation's coordinate
###################################################
#### Depends on different data ####
c = SkyCoord('11:26:57.771 +16:39:11.79', unit=(u.hourangle, u.deg))
###################################################
RA, DEC = c.ra.value, c.dec.value

def update_wcs(input, ext=0, dim=3): 
    iraf.hedit(input+'.fits['+str(ext)+']', 'CRVAL1', str(RA), add=yes, verify=no) 
    iraf.hedit(input+'.fits['+str(ext)+']', 'CRVAL2', str(DEC), add=yes, verify=no) 
    iraf.hedit(input+'.fits['+str(ext)+']', 'CRPIX1', str(CRPIX1), add=yes, verify=no)
    iraf.hedit(input+'.fits['+str(ext)+']', 'CRPIX2', str(CRPIX2), add=yes, verify=no)
    iraf.hedit(input+'.fits['+str(ext)+']', 'CD1_1', ' '+str( resample * np.cos((PA-180) / 180 * np.pi) ).replace('e-0','E-'), add=yes, verify=no) # dRA per 1 x-pixel
    iraf.hedit(input+'.fits['+str(ext)+']', 'CD2_1', ' '+str( resample * np.cos((PA-90) / 180 * np.pi) ).replace('e-0','E-'), add=yes, verify=no) # dDEC per 1 x-pixel
    iraf.hedit(input+'.fits['+str(ext)+']', 'CD1_2', ' '+str(resample * np.cos((PA-90) / 180 * np.pi) ).replace('e-0','E-'), add=yes, verify=no)
    iraf.hedit(input+'.fits['+str(ext)+']', 'CD2_2', ' '+str(resample * np.cos((PA) / 180 * np.pi) ).replace('e-0','E-'), add=yes, verify=no)
    iraf.hedit(input+'.fits['+str(ext)+']', 'CTYPE1', 'RA---TAN', add=yes, verify=no)
    iraf.hedit(input+'.fits['+str(ext)+']', 'CTYPE2', 'DEC--TAN', add=yes, verify=no)
    # double update the two items to avoid turncated decimals:
    iraf.hedit(input+'.fits['+str(ext)+']', 'CD2_1', ' '+str( resample * np.cos((PA-90) / 180 * np.pi) ).replace('e-0','E-'), add=yes, verify=no) # dDEC per 1 x-pixel
    iraf.hedit(input+'.fits['+str(ext)+']', 'CD1_2', ' '+str(resample * np.cos((PA-90) / 180 * np.pi) ).replace('e-0','E-'), add=yes, verify=no)
#     # convert in y-axis
#     ind = '[*,-*,*]' if dim==3 else '[*,-*]'
#     iraf.imcopy( input+'.fits['+str(ext)+']'+ind, input+'.fits['+str(ext)+', overwrite]' )
    # for showing in CARTA: (imcopy remove all CDELTs) 
    iraf.hedit(input+'.fits['+str(ext)+']', 'CDELT3', str(CD3_3), add=yes, verify=no)
    iraf.hedit(input+'.fits['+str(ext)+']', 'CTYPE1, CTYPE2, CRVAL1, CRVAL2, CRPIX1, CRPIX2, CD1_1, CD2_1, CD1_2, CD2_2, CDELT3', '.')

iraf.delete( 'wu'+prefix+obs+'.fits' )
iraf.copy( 'u'+prefix+obs+'.fits', 'wu'+prefix+obs+'.fits' )
for ext in range(5): update_wcs( 'wu'+prefix+obs, ext=ext, dim=3 )

os.system( 'gzip -f '+'wu'+prefix+obs+'.fits' ) # compress
# os.system( 'gzip -d '+'wu'+prefix+obs+'

# image
iraf.delete( 'w'+prefix+obs+'_wavestack.fits' )
iraf.copy( prefix+obs+'_wavestack.fits', 'w'+prefix+obs+'_wavestack.fits' )
update_wcs( 'w'+prefix+obs+'_wavestack', ext=0, dim=2 )

# delete tmp files
os.system( ' rm -v tmp* ' )

##############################################################
##############################################################

###### Determine PSF with acq images ###### 

iraf.psfmeasure( ','.join( [ '../observation/acq/N20210306S0047.fits['+str(id)+']' for id in [1,2,3,4] ] ), coords='markall', size='FWHM')
####
          Image  Column    Line     Mag    FWHM   Ellip      PA SAT
../observation/  180.03  577.03    1.45   4.188    0.03      47
                 170.66  922.20    2.62   4.252    0.01     -43
../observation/  224.47 1734.31    0.00   4.316    0.03      86
../observation/  112.28 1522.53    1.79   4.180    0.01     -87
                 161.73 1894.35    3.11   4.780    0.07     -72
../observation/  210.05 1676.10    2.45   4.228    0.12      18
#### follow J1126, 4.20*0.162=0.68

iraf.psfmeasure( ','.join( [ '../observation/acq/N20210322S0084.fits['+str(id)+']' for id in [1,2,3,4] ] ), coords='markall', size='FWHM')
####
          Image  Column    Line     Mag    FWHM   Ellip      PA SAT
../observation/   74.26 1470.74    6.38   4.552    0.09      -5
                 170.48  922.09    4.65   4.164    0.05      42
                 179.94  576.95    3.50   4.164    0.02      38
                 151.94  494.76    6.59   4.684    0.13      24
../observation/  100.80 1289.29    5.99   4.500    0.03       4
                 224.47 1734.03    2.04   4.364    0.01     -14
../observation/  161.96 1894.00    5.19   4.292    0.10     -62
                 112.32 1522.24    3.83   4.336    0.03      88
                 238.08  935.22    5.44   4.132    0.14     -32
                 177.98  786.90    6.87   2.826    0.14      78
../observation/  288.00 1421.89    4.03   2.252    0.37       0
                 210.15 1675.76    4.37   4.656    0.12      18
#### 4.20*0.162=0.68
