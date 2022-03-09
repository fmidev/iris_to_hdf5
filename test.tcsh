#!/bin/tcsh

# By this test script you can test the conversion from IRIS RAW data to ODIM HDF5.
# Just give the input filename and output filename, e.g.
# ./test.tcsh 201305272350_VAN.PPI2_A.raw test.h5

setenv ODIM_OUTPUT_FILE $2

setenv ODIM_OUTPUT_DIR .
setenv ODIM_COMPRESSION_LEVEL 6 # default 6, choose between 0 and 9
setenv ODIM_VOLUME_INTERVAL 5  # [min], nominal volume time is rounded using this 

setenv ODIM_Conventions 'ODIM_H5/V2_3'
setenv ODIM_what_version 'H5rad 2.3'
setenv ODIM_how_simulated False

setenv ODIM_ORIGCENTER 'EFKL' # WMO originating center

# Radar specific settings. The three letter radar site code is not included in the IRIS RAW
# metadata, but only in setup, so the code used here is the first three letters from 
# IRIS site name in RAW file. If not explicitely set, the code could be extracted from RAW file by command
# head --bytes=335 file.RAW | tail --bytes=3 



# Give quantities as ODIM-names (TH, DBZH, TV, DBZV, DBZHC, VRAD, VRADC, WRAD, SQI, 
# HCLASS, KDP, RHOHV, ZDR, LDR, TX, DBZX, SNR - see ODIM_struct.h)
# Adding number two to the end of name (like DBZH2) will produce 16-bit HDF5-data
# even if it's originally 8-bit. This works of course other way round also.

# Comma separated quantities could be asked, like 'DBZH,VRAD,HCLASS' 
# or '*' or ALL for all quantities present in RAW file. 

# Before the quantity name is the sweep number of IRIS volume, e.g. "3:DBZH" will
# convert DBZH data from the third sweep only. Comma-separated list gives several, and * all.

# You can combine whatever amount of IRIS RAW files to one HDF5 volume by converting 
# them first to *.dat files by IRIS_decoder and then give the list to the ODIM_encoder,
# e.g. ODIM_encoder F1.dat F2.data F3.dat 
# The output HDF5 file name must be given by setting the environment variable ODIM_OUTPUT_FILE.
# If the output file is not given, then the encoder makes ODIM-specified filename, e.g.
# T_PAGZ49_C_EFKL_20130531075002.h5

# You can add option -v as first argument to IRIS_decoder or ODIM_encoder
# to get list of conversions done (for e.g. log files)

setenv ODIM_OUTPUT_FILE $2

# Next variables can also be set per site, e.g. setenv ODIM_ESP_TXtype 'solid state'
setenv ODIM_TXtype 'magnetron'
setenv ODIM_poltype 'simultaneous-dual'
setenv ODIM_system 'VAISWRM200'

# Overall uptime reliability [%] if available:
# setenv ODIM_XXX_OUR 98

# Difference between calibration with continuous signal and pulsed operative mode
setenv FiniteBandwithLoss 1.4 
# From IRIS setup_dsp.conf:
# Antenna gain: ODIM_XXX_antgain = iantgain/10.0 

# IRIS_XXX_RXlossH = Horiz.iRcvLoss/10.0
# IRIS_XXX_RXlossV = Vert.iRcvLoss/10.0
# RXloss(ODIM) = RXloss(IRIS) - FiniteBandwithLoss
#                              
# When sending in H+V mode:
# TXlossH = Horiz.iXmtLossHAndV/10.0
# TXlossV = Vert.iXmtLossHAndV/10.0
# radconstHV(IRIS) = inghdr->tcf.cal.iRadarConstantVert/100.0
# (ingest_header -> task_configuration > task_calib_info -> iRadarConstantVert/100.0)
# radconstH(ODIM) = radconstHV(IRIS) + IRIS_XXX_RXlossH


setenv ODIM_KOR_source 'WIGOS:0-246-0-100926,WMO:02933,RAD:FI46,PLC:Korpo,NOD:fikor'
setenv IRIS_KOR_antgain 45.5
setenv IRIS_KOR_RXlossH 1.0
setenv IRIS_KOR_RXlossV 1.2
setenv IRIS_KOR_TXlossH 4.0
setenv IRIS_KOR_TXlossV 4.3
setenv ODIM_KOR_quantities '*:*'

setenv ODIM_VAN_source 'WIGOS:0-246-0-101001,WMO:02975,RAD:FI42,PLC:Vantaa,NOD:fivan'
setenv IRIS_VAN_antgain 45.2
setenv IRIS_VAN_RXlossH 1.0
setenv IRIS_VAN_RXlossV 1.2
setenv IRIS_VAN_TXlossH 4.0
setenv IRIS_VAN_TXlossV 4.3
setenv ODIM_VAN_quantities '*:*'

setenv ODIM_VIH_source 'WIGOS:0-246-0-000000,RAD:FI53,PLC:Vihti,NOD:fivih'
setenv IRIS_VIH_antgain 45.2
setenv IRIS_VIH_RXlossH 1.0
setenv IRIS_VIH_RXlossV 1.2
setenv IRIS_VIH_TXlossH 4.0
setenv IRIS_VIH_TXlossV 4.3
setenv ODIM_VIH_quantities '*:*'

setenv ODIM_ANJ_source 'WIGOS:0-246-0-101234,WMO:02954,RAD:FI44,PLC:Anjalankoski,NOD:fianj'
setenv IRIS_ANJ_antgain 45.8
setenv IRIS_ANJ_RXlossH 1.0
setenv IRIS_ANJ_RXlossV 1.2
setenv IRIS_ANJ_TXlossH 4.0
setenv IRIS_ANJ_TXlossV 4.3
setenv ODIM_ANJ_system 'VAISWRM200'
setenv ODIM_ANJ_quantities '*:*'

setenv ODIM_IKA_source 'WIGOS:0-246-0-101312,WMO:02942,RAD:FI43,PLC:Ikaalinen,NOD:fiika'
setenv IRIS_IKA_antgain 45.8
setenv IRIS_IKA_RXlossH 1.0
setenv IRIS_IKA_RXlossV 1.2
setenv IRIS_IKA_TXlossH 4.0
setenv IRIS_IKA_TXlossV 4.3
setenv ODIM_IKA_quantities '*:*'

setenv ODIM_KES_source 'WIGOS:0-246-0-100690,WMO:02995,RAD:FI50,PLC:Kesälahti,NOD:fikes'
setenv IRIS_KES_antgain 45.6
setenv IRIS_KES_RXlossH 1.0
setenv IRIS_KES_RXlossV 1.2
setenv IRIS_KES_TXlossH 4.0
setenv IRIS_KES_TXlossV 4.3
setenv ODIM_KES_quantities '*:*'

setenv ODIM_KUO_source 'WIGOS:0-246-0-101582,WMO:02918,RAD:FI45,PLC:Kuopio,NOD:fikuo'
setenv IRIS_KUO_antgain 45.0
setenv IRIS_KUO_RXlossH 1.0
setenv IRIS_KUO_RXlossV 1.2
setenv IRIS_KUO_TXlossH 4.0
setenv IRIS_KUO_TXlossV 4.3
setenv ODIM_KUO_quantities '*:*'

setenv ODIM_PET_source 'WIGOS:0-246-0-103813,WMO:02775,RAD:FI51,PLC:Petäjävesi,NOD:fipet'
setenv IRIS_PET_antgain 45.5
setenv IRIS_PET_RXlossH 1.0
setenv IRIS_PET_RXlossV 1.2
setenv IRIS_PET_TXlossH 4.0
setenv IRIS_PET_TXlossV 4.3
setenv ODIM_PET_quantities '*:*'

setenv ODIM_VIM_source 'WIGOS:0-246-0-101518,WMO:02925,RAD:FI49,PLC:Vimpeli,NOD:fivim'
setenv IRIS_VIM_antgain 45.1
setenv IRIS_VIM_RXlossH 1.0
setenv IRIS_VIM_RXlossV 1.2
setenv IRIS_VIM_TXlossH 4.0
setenv IRIS_VIM_TXlossV 4.3
setenv ODIM_VIM_quantities '*:*'

setenv ODIM_NUR_source 'WIGOS:0-246-0-107131,RAD:FI52,PLC:Nurmes,NOD:finur'
setenv IRIS_NUR_antgain 45.1
setenv IRIS_NUR_RXlossH 1.0
setenv IRIS_NUR_RXlossV 1.2
setenv IRIS_NUR_TXlossH 4.0
setenv IRIS_NUR_TXlossV 4.3
setenv ODIM_NUR_quantities '*:*'

setenv ODIM_UTA_source 'WIGOS:0-246-0-101872,WMO:02870,RAD:FI47,PLC:Utajärvi,NOD:fiuta'
setenv IRIS_UTA_antgain 45.1
setenv IRIS_UTA_RXlossH 1.0
setenv IRIS_UTA_RXlossV 1.2
setenv IRIS_UTA_TXlossH 4.0
setenv IRIS_UTA_TXlossV 4.3
setenv ODIM_UTA_quantities '*:*'

setenv ODIM_LUO_source 'WIGOS:0-246-0-101939,WMO:02840,RAD:FI48,PLC:Luosto,NOD:filuo'
setenv IRIS_LUO_antgain 47.5
setenv IRIS_LUO_RXlossH 1.0
setenv IRIS_LUO_RXlossV 1.2
setenv IRIS_LUO_TXlossH 4.0
setenv IRIS_LUO_TXlossV 4.3
setenv ODIM_LUO_quantities '*:*'

./IRIS_decoder $1 "$1".dat
./ODIM_encoder "$1".dat


