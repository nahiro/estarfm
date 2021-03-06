#!/usr/bin/env python
#---------------------------------------------------------------------------------
#                            ESTARFM PROGRAM
#               Using two pairs of fine and coarse images
#         the program can be used for whole TM scene and VI index product
#        Developed by (1) Zhu Xiaolin,email: zhuxiaolin55@gmail.com
#             Department of Land Surveying and Geo-Informatics
#             The Hong Kong Polytechnic University
#
#            Debugging history:
#            1)5/12/2012 correct the abnormal prediction
#            2)9/29/2013 correct the spatial distance calculation for integrate window
#            3)7/10/2014 correct the abnormal value of spectral distance and use all bands to indentify background
#            4)2/13/2017 add one parameter to specify the value of pixels of background or missing
#            5)1/1/2018  improve efficiency and modified the weight for fusing VI index
#            6)3/11/2008 correct a bug in spatial distance caculation
#            7)7/27/2018 correct a bug of abnormal coversion coefficent estimation when the two input pairs are too similar
#            8)7/27/2018 improve the prediction when no enough simular pixels are selected
#
# Please cite the reference: Xiaolin Zhu, Jin Chen, Feng Gao, & Jeffrey G Masek.
# An enhanced spatial and temporal adaptive reflectance fusion model for complex
# heterogeneous regions. Remote Sensing of Environment,2010,114,2610-2623
#
#                     Copyright belongs to Xiaolin Zhu
#---------------------------------------------------------------------------------

import os
import sys
from datetime import datetime
import gdal
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
from optparse import OptionParser,IndentedHelpFormatter

# ESTARFM

# please set the following parameters
#----------------------------------------------------------------------
HWID = 25.0           # set the haif window size, if 25, the window size is 25*2+1=51 fine pixels
NUM_CLASS = 4.0       # set the estimated number of classes, please set a larger value if blending images with very few bands
DN_MIN = 0            # set the range of DN value of the image, If byte, 0 and 255
DN_MAX = 10000.0      # set the range of DN value of the image, If byte, 0 and 255
BACKGROUND_F = -9999  # the value of background and missng pixels in fine images
BACKGROUND_C = -9999  # the value of background and missng pixels in coarse images
PATCH_LONG = 500      # set the size of each block,if process whole ETM scene, set 500-1000

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog fine_image_1 coarse_image_1 fine_image_2 coarse_image_2 coarse_image_0 [options]')
parser.add_option('-w','--hwid',default=HWID,type='float',help='set the haif window size, if 25, the window size is 25*2+1=51 fine pixels (%default)')
parser.add_option('-c','--num_class',default=NUM_CLASS,type='float',help='set the estimated number of classes, please set a larger value if blending images with very few bands (%default)')
parser.add_option('-m','--dn_min',default=DN_MIN,type='float',help='set the range of DN value of the image, If byte, 0 (%default)')
parser.add_option('-M','--dn_max',default=DN_MAX,type='float',help='set the range of DN value of the image, If byte, 255 (%default)')
parser.add_option('-b','--background_f',default=BACKGROUND_F,type='float',help='the value of background and missng pixels in fine images (%default)')
parser.add_option('-B','--background_c',default=BACKGROUND_C,type='float',help='the value of background and missng pixels in coarse images (%default)')
parser.add_option('-p','--patch_long',default=PATCH_LONG,type='int',help='set the size of each block,if process whole ETM scene, set 500-1000 (%default)')
parser.add_option('--fine_fname',default=None,help='File name of output fine image (%default)')
parser.add_option('--fine_dtype',default=None,help='Data type of output fine image, such as float32, uint16 (%default)')
(opts,args) = parser.parse_args()
if len(args) < 5:
    parser.print_help()
    sys.exit(0)
FileName1 = args[0]
FileName2 = args[1]
FileName3 = args[2]
FileName4 = args[3]
FileName5 = args[4]
if opts.fine_fname is None:
    opts.fine_fname = os.path.splitext(os.path.basename(FileName5))[0]+'_ESTARFM.tif'
wint = int(opts.hwid+0.1)

#-------------------------------------------------------------------
#                       main program
#-------------------------------------------------------------------

# open the fine image of the first pair
#------------------------------------------------------------------------
ds = gdal.Open(FileName1)
ns = ds.RasterXSize
nl = ds.RasterYSize
nb = ds.RasterCount
data_f1 = ds.ReadAsArray()
fine_trans = ds.GetGeoTransform()
fine_proj = ds.GetProjection()
ds = None
orig_ns = ns
orig_nl = nl
orig_nb = nb
if opts.fine_dtype is not None:
    fine_dtype = np.dtype(opts.fine_dtype)
else:
    fine_dtype = data_f1.dtype
n_ns = int(np.ceil(float(ns)/opts.patch_long)+0.1)
n_nl = int(np.ceil(float(nl)/opts.patch_long)+0.1)
n_sl = n_ns*n_nl
ind_patch = np.zeros((4,n_sl),dtype=np.intc)
for i_nl in range(n_nl):
    for i_ns in range(n_ns):
        ind_patch[0,n_ns*i_nl+i_ns] = i_ns*opts.patch_long
        ind_patch[1,n_ns*i_nl+i_ns] = min([ns,(i_ns+1)*opts.patch_long])
        ind_patch[2,n_ns*i_nl+i_ns] = i_nl*opts.patch_long
        ind_patch[3,n_ns*i_nl+i_ns] = min([nl,(i_nl+1)*opts.patch_long])
patch_f1 = []
for isub in range(n_sl):
    patch_f1.append(data_f1[:,ind_patch[2,isub]:ind_patch[3,isub],ind_patch[0,isub]:ind_patch[1,isub]])

# open the coarse image of the first pair
#-----------------------------------------------------------
ds = gdal.Open(FileName2)
if ds.RasterXSize != ns:
    raise ValueError('Error, RasterXSize={}, ns={} >>> {}'.format(ds.RasterXsize,ns,FileName2))
if ds.RasterYSize != nl:
    raise ValueError('Error, RasterYSize={}, nl={} >>> {}'.format(ds.RasterYsize,nl,FileName2))
if ds.RasterCount != nb:
    raise ValueError('Error, RasterCount={}, nb={} >>> {}'.format(ds.RasterCount,nb,FileName2))
data_c1 = ds.ReadAsArray()
ds = None
patch_c1 = []
for isub in range(n_sl):
    patch_c1.append(data_c1[:,ind_patch[2,isub]:ind_patch[3,isub],ind_patch[0,isub]:ind_patch[1,isub]])

# open the fine image of the second pair
#-----------------------------------------------------------
ds = gdal.Open(FileName3)
if ds.RasterXSize != ns:
    raise ValueError('Error, RasterXSize={}, ns={} >>> {}'.format(ds.RasterXsize,ns,FileName3))
if ds.RasterYSize != nl:
    raise ValueError('Error, RasterYSize={}, nl={} >>> {}'.format(ds.RasterYsize,nl,FileName3))
if ds.RasterCount != nb:
    raise ValueError('Error, RasterCount={}, nb={} >>> {}'.format(ds.RasterCount,nb,FileName3))
data_f2 = ds.ReadAsArray()
ds = None
patch_f2 = []
for isub in range(n_sl):
    patch_f2.append(data_f2[:,ind_patch[2,isub]:ind_patch[3,isub],ind_patch[0,isub]:ind_patch[1,isub]])

# open the coarse image of the second pair
#-----------------------------------------------------------
ds = gdal.Open(FileName4)
if ds.RasterXSize != ns:
    raise ValueError('Error, RasterXSize={}, ns={} >>> {}'.format(ds.RasterXsize,ns,FileName4))
if ds.RasterYSize != nl:
    raise ValueError('Error, RasterYSize={}, nl={} >>> {}'.format(ds.RasterYsize,nl,FileName4))
if ds.RasterCount != nb:
    raise ValueError('Error, RasterCount={}, nb={} >>> {}'.format(ds.RasterCount,nb,FileName4))
data_c2 = ds.ReadAsArray()
ds = None
patch_c2 = []
for isub in range(n_sl):
    patch_c2.append(data_c2[:,ind_patch[2,isub]:ind_patch[3,isub],ind_patch[0,isub]:ind_patch[1,isub]])

# open the coarse image of the prediction time
#-----------------------------------------------------------
ds = gdal.Open(FileName5)
if ds.RasterXSize != ns:
    raise ValueError('Error, RasterXSize={}, ns={} >>> {}'.format(ds.RasterXsize,ns,FileName5))
if ds.RasterYSize != nl:
    raise ValueError('Error, RasterYSize={}, nl={} >>> {}'.format(ds.RasterYsize,nl,FileName5))
if ds.RasterCount != nb:
    raise ValueError('Error, RasterCount={}, nb={} >>> {}'.format(ds.RasterCount,nb,FileName5))
data_c0 = ds.ReadAsArray()
ds = None
patch_c0 = []
for isub in range(n_sl):
    patch_c0.append(data_c0[:,ind_patch[2,isub]:ind_patch[3,isub],ind_patch[0,isub]:ind_patch[1,isub]])

#------------------------------------------------------------------
        #process  each block
#-------------------------------------------------------------------
t0 = datetime.now() # the initial time of program running

sys.stderr.write('there are total {} block(s)\n'.format(n_sl))

mosic_f0 = np.full((nb,nl,ns),opts.background_f,dtype=fine_dtype)

for isub in range(n_sl):

    # open each block image
    fine1 = patch_f1[isub]
    coarse1 = patch_c1[isub]
    fine2 = patch_f2[isub]
    coarse2 = patch_c2[isub]
    coarse0 = patch_c0[isub]

    nb,nl,ns = fine1.shape

    fine0 = np.zeros((nb,nl,ns)) # place the blended result

    # row index of images
    row_index = np.zeros((nl,ns),dtype=np.intc)
    for i in range(nl):
        row_index[i,:] = i
    # column index of images
    col_index = np.zeros((nl,ns),dtype=np.intc)
    for i in range(ns):
        col_index[:,i] = i

    # compute the uncertainty, 0.2% of each band is uncertain
    uncertain = (opts.dn_max*0.002)*(2**0.5)

    similar_th = np.zeros((2,nb)) # compute the threshold of similar pixel seeking

    for iband in range(nb):
        similar_th[0,iband] = np.nanstd(fine1[iband,:,:],ddof=1)*2.0/opts.num_class # pair 1
        similar_th[1,iband] = np.nanstd(fine2[iband,:,:],ddof=1)*2.0/opts.num_class # pair 2

    # compute the distance of each pixel in the window with the target pixel (integrate window)
    a = np.arange(opts.hwid*2.0+1.0).reshape(-1,1)
    b = np.ones_like(a).T
    D_D_all = 1.0+np.power(np.square(opts.hwid-np.dot(a,b))+np.square(opts.hwid-np.dot(b.T,a.T)),0.5)/opts.hwid

    # find interaction of valid pixels of all input images: exclude missing pixels and background
    cnd_valid = (fine1[0,:,:] != opts.background_f)
    cnd_valid &= (fine2[0,:,:] != opts.background_f)
    cnd_valid &= (coarse1[0,:,:] != opts.background_c)
    cnd_valid &= (coarse2[0,:,:] != opts.background_c)
    cnd_valid &= (coarse0[0,:,:] != opts.background_c)

    for j in range(nl): # retieve each target pixel

        for i in range(ns):

            if cnd_valid[j,i]: # do not process the background

                ai = max([0,i-wint]) # the window location
                bi = min([ns,i+wint+1])
                aj = max([0,j-wint])
                bj = min([nl,j+wint+1])

                cnd_wind_valid = cnd_valid[aj:bj,ai:bi]
                position_cand = np.full_like(cnd_wind_valid,True) # place the location of each similar pixel
                row_wind = row_index[aj:bj,ai:bi]
                col_wind = col_index[aj:bj,ai:bi]

                # searching for similar pixels
                for ipair in range(2):
                    for iband in range(nb):
                        if ipair == 0:
                            S_S = np.abs(fine1[iband,aj:bj,ai:bi]-fine1[iband,j,i])
                        else:
                            S_S = np.abs(fine2[iband,aj:bj,ai:bi]-fine2[iband,j,i])
                        position_cand &= (S_S < similar_th[ipair,iband])
                cnd_cand = position_cand & cnd_valid[aj:bj,ai:bi]
                number_cand = cnd_cand.sum()

                if (number_cand > 5):

                    x_cand = col_wind[cnd_cand]
                    y_cand = row_wind[cnd_cand]
                    finecand = np.zeros((number_cand,nb*2))
                    coasecand = np.zeros((number_cand,nb*2))
                    for ib in range(nb):
                        finecand[:,ib] = (fine1[ib,aj:bj,ai:bi])[cnd_cand]
                        finecand[:,ib+nb] = (fine2[ib,aj:bj,ai:bi])[cnd_cand]
                        coasecand[:,ib] = (coarse1[ib,aj:bj,ai:bi])[cnd_cand]
                        coasecand[:,ib+nb] = (coarse2[ib,aj:bj,ai:bi])[cnd_cand]

                    if (nb == 1): # for images with one band, like NDVI
                        S_D_cand = 1.0-0.5*(abs((finecand[:,0]-coasecand[:,0])/(finecand[:,0]+coasecand[:,0]))+abs((finecand[:,1]-coasecand[:,1])/(finecand[:,1]+coasecand[:,1]))) # compute the correlation
                    else:
                        # for images with multiple bands
                        sdx = np.nanstd(finecand,ddof=1,axis=1)
                        sdy = np.nanstd(coasecand,ddof=1,axis=1)
                        meanx = np.nanmean(finecand,axis=1)
                        meany = np.nanmean(coasecand,axis=1)
                        x_meanx = np.zeros((number_cand,nb*2))
                        y_meany = np.zeros((number_cand,nb*2))
                        for ib in range(nb*2):
                            x_meanx[:,ib] = finecand[:,ib]-meanx
                            y_meany[:,ib] = coasecand[:,ib]-meany
                        S_D_cand = nb*2.0*np.nanmean(x_meanx*y_meany,axis=1)/(sdx*sdy)/(nb*2.0-1)
                    cnd_nan = np.isnan(S_D_cand)
                    S_D_cand[cnd_nan] = 0.5 # correct the NaN value of correlation

                    if ((bi-ai)*(bj-aj) < (wint*2+1)*(wint*2+1)): # not an integrate window
                        D_D_cand = 1.0+np.power(np.square(i-x_cand)+np.square(j-y_cand),0.5)/opts.hwid # spatial distance
                    else:
                        D_D_cand = D_D_all[cnd_cand] # integrate window
                    C_D = (1.0-S_D_cand)*D_D_cand+0.0000001 # combined distance
                    weight = (1.0/C_D)/np.sum(1.0/C_D)

                    for iband in range(nb): # compute V
                        fine_cand = np.hstack(((fine1[iband,aj:bj,ai:bi])[cnd_cand],(fine2[iband,aj:bj,ai:bi])[cnd_cand]))
                        coarse_cand = np.hstack(((coarse1[iband,aj:bj,ai:bi])[cnd_cand],(coarse2[iband,aj:bj,ai:bi])[cnd_cand]))
                        coarse_change = abs(np.nanmean((coarse1[iband,aj:bj,ai:bi])[cnd_cand])-np.nanmean((coarse2[iband,aj:bj,ai:bi])[cnd_cand]))
                        if (coarse_change >= opts.dn_max*0.02): # to ensure changes in coarse image large enough to obtain the conversion coefficient
                            regress_result = sm.OLS(fine_cand,sm.add_constant(coarse_cand)).fit()
                            sig = 1.0-stats.f.cdf(regress_result.fvalue,1,number_cand*2-2)
                            # correct the result with no significancy or inconsistent change or too large value
                            if (sig <= 0.05) and (regress_result.params[1] > 0) and (regress_result.params[1] <= 5):
                                V_cand = regress_result.params[1]
                            else:
                                V_cand = 1.0
                        else:
                            V_cand = 1.0

                        # compute the temporal weight
                        difc_pair1 = np.abs(np.nanmean((coarse0[iband,aj:bj,ai:bi])[cnd_wind_valid])-np.nanmean((coarse1[iband,aj:bj,ai:bi])[cnd_wind_valid]))+0.01**5
                        difc_pair2 = np.abs(np.nanmean((coarse0[iband,aj:bj,ai:bi])[cnd_wind_valid])-np.nanmean((coarse2[iband,aj:bj,ai:bi])[cnd_wind_valid]))+0.01**5
                        T_weight1 = (1.0/difc_pair1)/(1.0/difc_pair1+1.0/difc_pair2)
                        T_weight2 = (1.0/difc_pair2)/(1.0/difc_pair1+1.0/difc_pair2)

                        # predict from pair1
                        coase0_cand = (coarse0[iband,aj:bj,ai:bi])[cnd_cand]
                        coase1_cand = (coarse1[iband,aj:bj,ai:bi])[cnd_cand]
                        fine01 = fine1[iband,j,i]+np.sum(weight*V_cand*(coase0_cand-coase1_cand))
                        # predict from pair2
                        coase2_cand = (coarse2[iband,aj:bj,ai:bi])[cnd_cand]
                        fine02 = fine2[iband,j,i]+np.sum(weight*V_cand*(coase0_cand-coase2_cand))
                        # the final prediction
                        fine0[iband,j,i] = T_weight1*fine01+T_weight2*fine02
                        # revise the abnormal prediction
                        if (fine0[iband,j,i] <= opts.dn_min) or (fine0[iband,j,i] >= opts.dn_max):
                            fine01 = np.sum(weight*(fine1[iband,aj:bj,ai:bi])[cnd_cand])
                            fine02 = np.sum(weight*(fine2[iband,aj:bj,ai:bi])[cnd_cand])
                            fine0[iband,j,i] = T_weight1*fine01+T_weight2*fine02
                else: # for the case of no enough similar pixel selected
                    for iband in range(nb):
                        # compute the temporal weight
                        difc_pair1 = np.nanmean((coarse0[iband,aj:bj,ai:bi])[cnd_wind_valid])-np.nanmean((coarse1[iband,aj:bj,ai:bi])[cnd_wind_valid])+0.01**5
                        difc_pair1_a = np.abs(difc_pair1)
                        difc_pair2 = np.nanmean((coarse0[iband,aj:bj,ai:bi])[cnd_wind_valid])-np.nanmean((coarse2[iband,aj:bj,ai:bi])[cnd_wind_valid])+0.01**5
                        difc_pair2_a = np.abs(difc_pair2)
                        T_weight1 = (1.0/difc_pair1_a)/(1.0/difc_pair1_a+1.0/difc_pair2_a)
                        T_weight2 = (1.0/difc_pair2_a)/(1.0/difc_pair1_a+1.0/difc_pair2_a)
                        fine0[iband,j,i] = T_weight1*(fine1[iband,j,i]+difc_pair1)+T_weight2*(fine2[iband,j,i]+difc_pair2)

    # mosiac the blended patch
    mosic_f0[:,ind_patch[2,isub]:ind_patch[3,isub],ind_patch[0,isub]:ind_patch[1,isub]] = fine0.astype(fine_dtype) # change the type of prediction into the type same as the input image
    sys.stderr.write('finished {} block\n'.format(isub+1))

#--------------------------------------------------------------------------------------
# Output result
if fine_dtype is np.dtype('int16'):
    output_dtype = gdal.GDT_Int16
elif fine_dtype is np.dtype('uint16'):
    output_dtype = gdal.GDT_UInt16
elif fine_dtype is np.dtype('int32'):
    output_dtype = gdal.GDT_Int32
elif fine_dtype is np.dtype('uint32'):
    output_dtype = gdal.GDT_UInt32
elif fine_dtype is np.dtype('float32'):
    output_dtype = gdal.GDT_Float32
elif fine_dtype is np.dtype('float64'):
    output_dtype = gdal.GDT_Float64
else:
    raise ValueError('Error in data type >>> {}'.format(fine_dtype))
driver = gdal.GetDriverByName("GTiff")
ds = driver.Create(opts.fine_fname,orig_ns,orig_nl,orig_nb,output_dtype)
ds.SetGeoTransform(fine_trans) # sets same geotransform as input
ds.SetProjection(fine_proj) # sets same projection as input
for iband in range(nb):
    band = ds.GetRasterBand(iband+1)
    band.WriteArray(mosic_f0[iband])
band.SetNoDataValue(np.nan) # if you want these values transparent
ds.FlushCache() # saves to disk
ds = None

t1 = datetime.now()
dt = (t1-t0).total_seconds()
sys.stderr.write('time used: {:.0f} h {:.0f} m {:.1f} s\n'.format(dt/3600,np.mod(dt,3600)/60,np.mod(dt,60)))
