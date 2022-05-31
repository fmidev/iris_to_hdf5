#!/bin/bash

# By this test script you can test the conversion from IRIS RAW data to ODIM HDF5.
# Just give the input filename and output filename, e.g.
# ./test.bash 201305272350_VAN.PPI2_A.raw test.h5
echo -e "\n#########################################################################"
echo -e "\nRunning conversion tests for IRIS_decoder and ODIM_encoder\n\n"

DECODER=bin/IRIS_decoder
ENCODER=bin/ODIM_encoder

RAW=testdata/201303151250_VAN.PPI2_E.raw

export ODIM_OUTPUT_FILE=test.h5
export ODIM_OUTPUT_DIR=.
export ODIM_COMPRESSION_LEVEL=6 # default 6, choose between 0 and 9
export ODIM_VOLUME_INTERVAL=5  # [min], nominal volume time is rounded using this 

export ODIM_Conventions='ODIM_H5/V2_3'
export ODIM_what_version='H5rad 2.3'
export ODIM_how_simulated=False

export ODIM_ORIGCENTER='EFKL' # WMO originating center
# export ODIM_ORIGCENTER='EEMH' # WMO originating center

# Radar specific settings. The three letter radar site code is not included in the IRIS RAW
# metadata, but only in setup, so the code used here is the first three letters from 
# IRIS site name in RAW file.


# Next variables can also be set per site, e.g. export ODIM_ESP_TXtype 'solid state'
export ODIM_VAN_source='WIGOS:0-246-0-101001,WMO:02975,RAD:FI42,PLC:Vantaa,NOD:fivan'
export ODIM_VAN_system='VAISWRM200'
export ODIM_TXtype='magnetron'
export ODIM_poltype='simultaneous-dual'
export ODIM_system='VAISWRM200'

# Overall uptime reliability [%] if available:
# export ODIM_XXX_OUR 98

# Difference between calibration with continuous signal and pulsed operative mode
export FiniteBandwithLoss=1.4 
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

export IRIS_VAN_antgain=45.2
export IRIS_VAN_RXlossH=1.0
export IRIS_VAN_RXlossV=1.2
export IRIS_VAN_TXlossH=4.0
export IRIS_VAN_TXlossV=4.3

# Give quantities as ODIM-names (TH, DBZH, TV, DBZV DBZHC, VRAD, VRADC, WRAD, SQI, 
# HCLASS, KDP, RHOHV, ZDR, LDR, TX, DBZX - see ODIM_struct.h)
# Adding number two to the end of name (like DBZH2) will produce 16-bit HDF5-data
# even if it's originally 8-bit. This works of course other way round also.

# Comma-separated several quantities could be asked, like 'DBZH,VRAD,HCLASS' 
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
# The -d option of IRIS_decoder dumps all metadata

export ODIM_VAN_quantities='*:*'

${DECODER} -v  $RAW "$RAW".dat
DECRES=$?
${ENCODER} -v "$RAW".dat
ENCRES=$?
rm "$RAW".dat
ODIM_PATH=$ODIM_OUTPUT_DIR/$ODIM_OUTPUT_FILE
echo -e "\nIRIS RAW $RAW conversion to $ODIM_PATH\n${DECODER} return ${DECRES}\n${ENCODER} return ${ENCRES}\n"

exit
