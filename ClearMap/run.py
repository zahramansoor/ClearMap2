#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 11:52:21 2021

@author: kepecs

zmd's makeshift pipeline mimicking ClearMap2/CellMap
from jupyter tutorial
"""

#ClearMap path
import sys
sys.path.append('/home/kepecs/ClearMap2')

#load ClearMap modules
from ClearMap.Environment import *  #analysis:ignore

#directories and files
directory = "/media/kepecs/Starosta_8T/2P_imaging/AA6-AK1a"    

expression_raw      = 'AA6-AK1a_640/AA6-AK1a_640-stitched_T001_Z<Z,I>_C01.tif'            
expression_auto     = 'AA6-AK1a_561/original/AA6-AK1a_561-stitched_T001_Z<Z,I>_C01.tif'

ws = wsp.Workspace('CellMap', directory=directory);
ws.update(raw=expression_raw, autofluorescence=expression_auto)
ws.debug = False

resources_directory = settings.resources_path

ws.info()

#init atals and reference files
annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
    slicing=(slice(None),slice(None),slice(0,256)), orientation=(1,-2,3),
    overwrite=False, verbose=True);

#alignment parameter files    
align_channels_affine_file   = io.join(resources_directory, 'Alignment/align_affine.txt')
align_reference_affine_file  = io.join(resources_directory, 'Alignment/align_affine.txt')
align_reference_bspline_file = io.join(resources_directory, 'Alignment/align_bspline.txt')

#%% Convet raw data to MEMORY MAPPED npy file     
             
source = ws.source('autofluorescence');
sink   = ws.filename('stitched', postfix = 'autofluorescence')
io.convert(source, sink, verbose=True)

resample_parameter = {
    "source_resolution" : (1.625,1.625,1.6),
    "sink_resolution"   : (25,25,25),
    "orientation": (-3, -2, 1), #inverts old z (dorsal -> ventral) and y and flips x and z
    "processes" : None,
    "verbose" : True,             
    };

res.resample(ws.filename('stitched', postfix = 'autofluorescence'), sink=ws.filename('resampled', postfix = 'autofluorescence'), **resample_parameter)