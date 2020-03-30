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

from datetime import datetime
import gdal
import numpy as np
from optparse import OptionParser,IndentedHelpFormatter

# ESTARFM

# please set the following parameters
#----------------------------------------------------------------------
hwid = 25.0              # set the haif window size, if 25, the window size is 25*2+1=51 fine pixels
num_class = 4.0       # set the estimated number of classes, please set a larger value if blending images with very few bands
dn_min = 0            # set the range of DN value of the image,If byte, 0 and 255
dn_max = 10000.0
background_f = -9999  # the value of background and missng pixels in fine images
background_c = -9999  # the value of background and missng pixels in coarse images
patch_long = 500      # set the size of each block,if process whole ETM scene, set 500-1000
temp_file = '/tmp'    # set the temporary file location, temporary files will be deleted after the work

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog fine_image_1 coarse_image_1 fine_image_2 coarse_image_2 coarse_image_0 [options]')
(opts,args) = parser.parse_args()
if len(args) < 5:
    parser.print_help()
    sys.exit(0)
FileName1 = args[0]
FileName2 = args[1]
FileName3 = args[2]
FileName4 = args[3]
FileName5 = args[4]
wint = int(hwid+0.1)

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
ds = None
orig_ns = ns
orig_nl = nl
n_ns = int(np.ceil(float(ns)/patch_long)+0.1)
n_nl = int(np.ceil(float(nl)/patch_long)+0.1)
n_sl = n_ns*n_nl
ind_patch = np.zeros((4,n_sl))
for i_nl in range(n_nl):
    for i_ns in range(n_ns):
        ind_patch[0,n_ns*i_nl+i_ns] = i_ns*patch_long
        ind_patch[1,n_ns*i_nl+i_ns] = min([ns,(i_ns+1)*patch_long])
        ind_patch[2,n_ns*i_nl+i_ns] = i_nl*patch_long
        ind_patch[3,n_ns*i_nl+i_ns] = min([nl,(i_nl+1)*patch_long])
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

print('there are total',n_sl,' blocks')

for isub in range(n_sl):

    # open each block image
    fine1 = data_f1[isub]
    coarse1 = data_c1[isub]
    fine2 = data_f2[isub]
    coarse2 = data_c2[isub]
    coarse0 = data_c0[isub]

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
    uncertain = (dn_max*0.002)*(2**0.5)

    similar_th = fltarr(2,nb) # compute the threshold of similar pixel seeking

    for iband in range(nb):
        similar_th[0,iband] = np.nanstd(fine1[iband,:,:])*2.0/num_class #pair 1
        similar_th[1,iband] = np.nanstd(fine2[iband,:,:])*2.0/num_class #pair 2

    # compute the distance of each pixel in the window with the target pixel (integrate window)
    a = np.arange(hwid*2.0+1.0).reshape(-1,1)
    b = np.ones_like(a).T
    D_D_all = 1.0+np.power(np.square(hwid-np.dot(a,b))+np.square(hwid-np.dot(b.T,a.T)),0.5)/hwid
    
    # find interaction of valid pixels of all input images: exclude missing pixels and background
    cnd_valid = (fine1[0,:,:] != background_f)
    cnd_valid &= (fine2[0,:,:] != background_f)
    cnd_valid &= (coarse1[0,:,:] != background_c)
    cnd_valid &= (coarse2[0,:,:] != background_c)
    cnd_valid &= (coarse0[0,:,:] != background_c)
     
    for j in range(nl): # retieve each target pixel

        for i in range(ns):
       
            if cnd_valid[j,i]: # do not process the background

                ai = max([0,i-wint]) # the window location
                bi = min([ns,i+wint])
                aj = max([0,j-wint])
                bj = min([nl,j+wint])

                cnd_wind_valid = cnd_valid[aj:bj,ai:bi]
                position_cand = np.full((bi-ai+1)*(bj-aj+1),True) # place the location of each similar pixel
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

                    S_D_cand = fltarr(number_cand)                #compute the correlation
                    x_cand = col_wind[cnd_cand]
                    y_cand = row_wind[cnd_cand]
                    finecand = fltarr(number_cand,nb*2)
                    coasecand = fltarr(number_cand,nb*2)
                    for ib=0,nb-1, 1 do begin
                        finecand[*,ib] = (fine1[ib,aj:bj,ai:bi])[cnd_cand]
                        finecand[*,ib+nb] = (fine2[ib,aj:bj,ai:bi])[cnd_cand]
                        coasecand[*,ib] = (coarse1[ib,aj:bj,ai:bi])[cnd_cand]
                        coasecand[*,ib+nb] = (coarse2[ib,aj:bj,ai:bi])[cnd_cand]
                    endfor
                 
                    if (nb eq 1) then begin  # for images with one band, like NDVI
                        S_D_cand=1.0-0.5*(abs((finecand[*,0]-coasecand[*,0])/(finecand[*,0]+coasecand[*,0]))+abs((finecand[*,1]-coasecand[*,1])/(finecand[*,1]+coasecand[*,1])))
                    endif else begin   
                        # for images with multiple bands             
                        sdx=stddev(finecand,DIMENSION=2)
                        sdy=stddev(coasecand,DIMENSION=2)         
                        meanx=mean(finecand,DIMENSION=2)
                        meany=mean(coasecand,DIMENSION=2)
                        x_meanx=fltarr(number_cand,nb*2)
                        y_meany=fltarr(number_cand,nb*2)
                        for ib=0,nb*2-1, 1 do begin
                            x_meanx[*,ib]=finecand[*,ib]-meanx
                            y_meany[*,ib]=coasecand[*,ib]-meany
                        endfor     
                        S_D_cand=nb*2.0*mean(x_meanx*y_meany,DIMENSION=2)/(sdx*sdy)/(nb*2.0-1) 
                    endelse
                    ind_nan = where(S_D_cand ne S_D_cand,num_nan)
                    if (num_nan gt 0) then S_D_cand[ind_nan]=0.5 # correct the NaN value of correlation

                    D_D_cand = fltarr(number_cand) # spatial distance
                    if ((bi-ai+1)*(bj-aj+1) lt (w*2.0+1)*(w*2.0+1)) then begin # not an integrate window
                         D_D_cand = 1.0+((i-x_cand)^2+(j-y_cand)^2)^0.5/float(w)              
                    endif else begin
                         D_D_cand[0:number_cand-1] = D_D_all[cnd_cand] # integrate window
                    endelse
                    C_D = (1.0-S_D_cand)*D_D_cand+0.0000001 # combined distance
                    weight = (1.0/C_D)/total(1.0/C_D)

                    for iband=0,nb-1,1 do begin # compute V
                        fine_cand=[(fine1[iband,aj:bj,ai:bi])[cnd_cand],(fine2[iband,aj:bj,ai:bi])[cnd_cand]]
                        corse_cand=[(coarse1[iband,aj:bj,ai:bi])[cnd_cand],(coarse2[iband,aj:bj,ai:bi])[cnd_cand]]
                        coarse_change=abs(mean((coarse1[iband,aj:bj,ai:bi])[cnd_cand])-mean((coarse2[iband,aj:bj,ai:bi])[cnd_cand]))         
                        if ( coarse_change ge dn_max*0.02) then begin #to ensure changes in coarse image large enough to obtain the conversion coefficient
                            regress_result=regress(corse_cand,fine_cand,FTEST=fvalue)
                            sig = 1.0-f_pdf(fvalue,1,number_cand*2-2)
                            # correct the result with no significancy or inconsistent change or too large value  
                            if (sig le 0.05 and regress_result[0] gt 0 and regress_result[0] le 5) then begin
                                V_cand = regress_result[0]
                            endif else begin
                                V_cand = 1.0
                            endelse
                        endif else begin
                            V_cand = 1.0
                        endelse

                        # compute the temporal weight
                        difc_pair1 = abs(mean((coarse0[iband,aj:bj,ai:bi])[ind_wind_valid])-mean((coarse1[iband,aj:bj,ai:bi])[ind_wind_valid]))+0.01^5
                        difc_pair2 = abs(mean((coarse0[iband,aj:bj,ai:bi])[ind_wind_valid])-mean((coarse2[iband,aj:bj,ai:bi])[ind_wind_valid]))+0.01^5
                        T_weight1 = (1.0/difc_pair1)/(1.0/difc_pair1+1.0/difc_pair2)
                        T_weight2 = (1.0/difc_pair2)/(1.0/difc_pair1+1.0/difc_pair2)

                        # predict from pair1
                        coase0_cand = (coarse0[iband,aj:bj,ai:bi])[cnd_cand]
                        coase1_cand = (coarse1[iband,aj:bj,ai:bi])[cnd_cand]
                        fine01 = fine1[iband,j,i]+total(weight*V_cand*(coase0_cand-coase1_cand))
                        # predict from pair2
                        coase2_cand = (coarse2[iband,aj:bj,ai:bi])[cnd_cand]
                        fine02 = fine2[iband,j,i]+total(weight*V_cand*(coase0_cand-coase2_cand))
                        # the final prediction
                        fine0[iband,j,i] = T_weight1*fine01+T_weight2*fine02
                        # revise the abnormal prediction
                        if (fine0[iband,j,i] le dn_min or fine0[iband,j,i] ge dn_max) then begin
                                fine01 = total(weight*(fine1[iband,aj:bj,ai:bi])[cnd_cand])
                                fine02 = total(weight*(fine2[iband,aj:bj,ai:bi])[cnd_cand])  
                                fine0[iband,j,i] = T_weight1*fine01+T_weight2*fine02
                        endif
                    endfor
                endif else begin   # for the case of no enough similar pixel selected
                    for iband=0,nb-1,1 do begin  
                        # compute the temporal weight
                        difc_pair1 = mean((coarse0[iband,aj:bj,ai:bi])[ind_wind_valid])-mean((coarse1[iband,aj:bj,ai:bi])[ind_wind_valid])+0.01^5
                        difc_pair1_a = abs(difc_pair1)
                        difc_pair2 = mean((coarse0[iband,aj:bj,ai:bi])[ind_wind_valid])-mean((coarse2[iband,aj:bj,ai:bi])[ind_wind_valid])+0.01^5
                        difc_pair2_a = abs(difc_pair2)
                        T_weight1 = (1.0/difc_pair1_a)/(1.0/difc_pair1_a+1.0/difc_pair2_a)
                        T_weight2 = (1.0/difc_pair2_a)/(1.0/difc_pair1_a+1.0/difc_pair2_a)
                        fine0[iband,j,i] = T_weight1*(fine1[iband,j,i]+difc_pair1)+T_weight2*(fine2[iband,j,i]+difc_pair2)
                    endfor
                endelse
            endif
        endfor
    endfor

    # change the type of prediction into the type same as the input image
    case Data_Type Of
        1:fine0 = Byte(fine0)    #  BYTE  Byte
        2:fine0 = FIX(fine0)     #  INT  Integer
        3:fine0 = LONG(fine0)    #  LONG  Longword integer
        4:fine0 = FLOAT(fine0)   #  FLOAT  Floating point
        5:fine0 = DOUBLE(fine0)  #  DOUBLE  Double-precision floating
        6:fine0 = COMPLEX(fine0)# complex, single-precision, floating-point
        9:fine0 = DCOMPLEX(fine0)#complex, double-precision, floating-point
        12:fine0 = UINT(fine0)   # unsigned integer vector or array
        13:fine0 = ULONG(fine0)   #  unsigned longword integer vector or array
        14:fine0 = LONG64(fine0)   #a 64-bit integer vector or array
        15:fine0a = ULONG64(fine0)   #an unsigned 64-bit integer vector or array
    EndCase

    print,'finished ',isub+1,' block'
    tempoutname1=temp_file+'/temp_blended'
    Envi_Write_Envi_File,fine0,Out_Name = tempoutname1+strtrim(isub+1,1)
    envi_file_mng, id=Fid1, /remove, /delete
    envi_file_mng, id=Fid2, /remove, /delete
    envi_file_mng, id=Fid3, /remove, /delete
    envi_file_mng, id=Fid4, /remove, /delete
    envi_file_mng, id=Fid5, /remove, /delete
endfor

#--------------------------------------------------------------------------------------
# mosiac all the blended patch

  mfid=intarr(n_sl)
  mdims=intarr(5,n_sl)
  mpos=intarr(nb,n_sl)
  pos=indgen(nb)
  x0=intarr(n_sl)
  y0=intarr(n_sl)

  for isub=0,n_sl-1,1 do begin
      envi_open_file, tempoutname1+strtrim(isub+1,1), r_fid= sub_fid
     if (sub_fid eq -1) then begin
       envi_batch_exit
       return
     endif
      envi_file_query,  sub_fid, ns=sub_ns, nl=sub_nl
      mfid[isub] = sub_fid
      mpos[*,isub] = indgen(nb)
      mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
      x0[isub]=ind_patch[0,isub]
      y0[isub]=ind_patch[2,isub]
  endfor

    xsize = orig_ns
    ysize = orig_nl
    pixel_size = [1.,1.]

    use_see_through = replicate(1L,n_sl)
    see_through_val = replicate(0L,n_sl)

    out_name=FileName5+'_ESTARFM'
    envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
    dims=mdims, out_name=out_name, xsize=xsize, $
    ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
    out_dt=Data_Type, pixel_size=pixel_size, $
    background=0, see_through_val=see_through_val, $
    use_see_through=use_see_through

    for i=0,n_sl-1,1 do begin
      envi_file_mng, id=mfid[i], /remove, /delete
    endfor

t1 = datetime.now()
dt = (t1-t0).totalsecond()
print('time used:', dt//3600,'h',mod(dt,3600)//60,'m',mod(dt,60),'s')
