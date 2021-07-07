#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CellMap
=======

This script is the main pipeline to analyze immediate early gene expression 
data from iDISCO+ cleared tissue [Renier2016]_.

See the :ref:`CellMap tutorial </CellMap.ipynb>` for a tutorial and usage.


.. image:: ../Static/cell_abstract_2016.jpg
   :target: https://doi.org/10.1016/j.cell.2020.01.028
   :width: 300

.. figure:: ../Static/CellMap_pipeline.png

  iDISCO+ and ClearMap: A Pipeline for Cell Detection, Registration, and 
  Mapping in Intact Samples Using Light Sheet Microscopy.


References
----------
.. [Renier2016] `Mapping of brain activity by automated volume analysis of immediate early genes. Renier* N, Adams* EL, Kirst* C, Wu* Z, et al. Cell. 2016 165(7):1789-802 <https://doi.org/10.1016/j.cell.2016.05.007>`_
"""
__author__    = 'Christoph Kirst <christoph.kirst.ck@gmail.com>'
__license__   = 'GPLv3 - GNU General Pulic License v3 (see LICENSE)'
__copyright__ = 'Copyright © 2020 by Christoph Kirst'
__webpage__   = 'http://idisco.info'
__download__  = 'http://www.github.com/ChristophKirst/ClearMap2'

#ClearMap path
import sys, tifffile as tif, os
import argparse
sys.path.append('/home/kepecs/python/ClearMap2/')

def parsecmdline():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "brainname", help="Brain name")
    parser.add_argument(
        "zstep", help="Z step size (in um)",
        type=float)
    parser.add_argument(
        "--channel1", help="Cell channel 1", type=int, default=640)
    parser.add_argument(
        "--channel2", help="Cell channel 2", type=int, default=561) #allows for cell detection on more than 1 channel per brain
    parser.add_argument(
        "--backgroundparam", help="Background subtraction parameter",
        type=tuple, default=(10,10))
    parser.add_argument(
        "--cellshape", help="Cell shape parameter",
        type=int, default=500)
    parser.add_argument(
        "--maximashape", help="Maxima detection shape",
        type=int, default=5)
    parser.add_argument(
        "--orientation", help="Reorientation paramters in xyz, e.g. (-3,-2,1) flips z and y",
        type=tuple, default=(-3, -2, 1))
    parser.add_argument(
        "--dst", help="Save destination",
        type=str, default="/home/kepecs/Documents/")
    return parser.parse_args()

def fillargs(args):
    """kepecs lab experiment specific params"""
    args.src = "/mnt/uncertainty"
    args.background640 = (10,10)
    args.cellshape640 = 500
    args.maximashape640 = 5
    args.background561 = None
    args.cellshape561 = 2600
    args.maximashape561 = 10
    
    return args

if __name__ == "__main__":
     
  #%############################################################################
  ### Initialization 
  ##############################################################################
  
  #% Initialize workspace
  from ClearMap.Environment import *  #analysis:ignore
  
  #directories and files
  args = parsecmdline()
  args = fillargs(args)
  directory = os.path.join(args.src, args.brainname)
  #identify strings to feed into regex
  str_auto = os.listdir(os.path.join(directory, "{0}_{1}".format(args.brainname, 488)))[98][:-12]
  str_ch1 = os.listdir(os.path.join(directory, "{0}_{1}".format(args.brainname, args.channel1)))[89][:-12]
  expression_ch1      = "{0}_{1}/{2}Z<Z,I>_C01.tif".format(args.brainname, args.channel1, str_ch1)
  expression_auto     = "{0}_{1}/{2}Z<Z,I>_C01.tif".format(args.brainname, 488, str_auto)

  ws = wsp.Workspace("CellMap", directory=directory);
  if isinstance(args.channel2,int):
      str_ch2 = os.listdir(os.path.join(directory, "{0}_{1}".format(args.brainname, args.channel2)))[58][:-12]
      expression_ch2      = "{0}_{1}/{2}Z<Z,I>_C01.tif".format(args.brainname, args.channel2, str_ch2)
      ws.update(ch1=expression_ch1, ch2=expression_ch2, autofluorescence=expression_auto)
  else:
      ws.update(ch1=expression_ch1, autofluorescence=expression_auto)
  ws.debug = False
    
  resources_directory = settings.resources_path
    
  ws.info()
      
  resources_directory = settings.resources_path
  
  #%% Initialize alignment 

  #alignment parameter files    
  align_channels_affine_file   = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_affine_file  = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_bspline_file = io.join(resources_directory, 'Alignment/Order2_Par0000bspline.txt')
  
  #%%############################################################################
  ### Resampling and atlas alignment 
  ###############################################################################
        
  #%% Resample 
  
  resample_parameter = {
        "source_resolution" : (1.26, 1.26, args.zstep), #z step might be 5 because less z planes
        "sink_resolution"   : (25,25,25),
        "orientation": args.orientation, #inverts old z (dorsal -> ventral) and flips x and z
        "processes" : None,
        "verbose" : True,             
        };
  #mod because of server priviledge issues
  sink_ch1 = os.path.join(args.dst, args.brainname+"_{}_resampled.tif".format(args.channel1))   
  # res.resample(ws.source("ch1"), sink=sink_ch1, **resample_parameter) 
  if isinstance(args.channel2,int):
      sink_ch2 = os.path.join(args.dst, args.brainname+"_{}_resampled.tif".format(args.channel2))   
      # res.resample(ws.source("ch2"), sink=sink_ch2, **resample_parameter) 

  #%% Resample autofluorescence
      
  resample_parameter_auto = {
        "source_resolution" : (1.26, 1.26, args.zstep), #z step might be 5 because less z planes
        "sink_resolution"   : (25,25,25),
        "orientation": args.orientation, #inverts old z (dorsal -> ventral) and flips x and z
        "processes" : None,
        "verbose" : True,             
        };    
  #mod because of server priviledge issues
  sink_auto = os.path.join(args.dst, args.brainname+"_auto_resampled.tif")   
  # res.resample(ws.filename("autofluorescence"), sink=sink_auto, **resample_parameter_auto) 
  # ws.update(resampled = sink_auto, resampled_to_auto_ch1 = sink_ch1)
  # if isinstance(args.channel2,int):
      # ws.update(resampled_to_auto_ch2 = sink_ch2)
  #%% Aignment - resampled to autofluorescence
  
  #align the two channels
  align_channels_parameter = {            
       #moving and reference images
       "moving_image" : sink_auto,
       "fixed_image"  : sink_ch1,
      
       #elastix parameter files for alignment
       "affine_parameter_file"  : align_channels_affine_file,
       "bspline_parameter_file" : None,
      
       "result_directory" :  os.path.join(args.dst, args.brainname+"_elastix_resampled_to_auto_ch1" )
       }; 
  
  # elx.align(**align_channels_parameter);
  
  if isinstance(args.channel2,int):
      # align the two channels
      align_channels_parameter = {            
          #moving and reference images
          "moving_image" : sink_auto,
          "fixed_image"  : sink_ch2,
          
          #elastix parameter files for alignment
          "affine_parameter_file"  : align_channels_affine_file,
          "bspline_parameter_file" : None,
          
          "result_directory" :  os.path.join(args.dst, args.brainname+"_elastix_resampled_to_auto_ch2" )
          }; 
      
      # elx.align(**align_channels_parameter);
  
  
  #%% Alignment - autoflourescence to reference
  
  #align autofluorescence to reference
  align_reference_parameter = {            
       #moving and reference images
       "moving_image" : "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_reference.tif", #whole brain
       "fixed_image"  : sink_auto,
      
       #elastix parameter files for alignment
       "affine_parameter_file"  :  align_reference_affine_file, 
       "bspline_parameter_file" :  align_reference_bspline_file, #mods from brainpipe
       #directory of the alignment result
       "result_directory" :  os.path.join(args.dst, args.brainname+"_elastix_auto_to_reference")
       };
  
  # elx.align(**align_reference_parameter);
  

  #%%############################################################################
  ### Cell detection
  ###############################################################################
  
  #%% Cell detection:
  
  #channel 1
  if args.channel1 == 640:
      background = args.background640
      cellshape = args.cellshape640
      maximashape = args.maximashape640
  elif args.channel1 == 561:
      background = args.background561
      cellshape = args.cellshape561
      maximashape = args.maximashape561
      
  cell_detection_parameter = cells.default_cell_detection_parameter.copy();
  cell_detection_parameter["illumination_correction"] = None;
  if background is not None: 
      cell_detection_parameter["background_correction"] = {"shape": background, "form": "Disk"}; 
  else: cell_detection_parameter["background_correction"] = background
  cell_detection_parameter["intensity_detection"]["measure"] = ["source"];
  cell_detection_parameter["shape_detection"]["threshold"] = cellshape #for 640 ch 500
  cell_detection_parameter["maxima_detection"]["shape"] = maximashape
  cell_detection_parameter["maxima_detection"]["save"] = False #DO NOT SAVE MAXIMA WTF

  processing_parameter = cells.default_cell_detection_processing_parameter.copy();
  processing_parameter.update(
        processes = 1, # 'serial', #multiple processes don't work on kepecs desktop bc of memory
        size_max = 10, #100, #35,
        size_min = 5,# 30, #30,
        optimization = False,
        optimization_fix = None,
        overlap  = 3, #32, #10,
        verbose = True
        )
  dstch1 = os.path.join(args.dst, "{0}_cells_{1}_raw.npy".format(args.brainname, args.channel1)) 
  print("\n\n STARTING CELL DETECTION FOR CHANNEL: {0}".format(args.channel1))   
  cells.detect_cells(ws.source("ch1"), dstch1,
                       cell_detection_parameter=cell_detection_parameter, 
                       processing_parameter=processing_parameter)
  
  #channel 2
  if args.channel2 == 640:
      background = args.background640
      cellshape = args.cellshape640
      maximashape = args.maximashape640
  elif args.channel2 == 561:
      background = args.background561
      cellshape = args.cellshape561
      maximashape = args.maximashape561
      
  cell_detection_parameter = cells.default_cell_detection_parameter.copy();
  cell_detection_parameter["illumination_correction"] = None;
  if background is not None: 
      cell_detection_parameter["background_correction"] = {"shape": background, "form": "Disk"}; 
  else: cell_detection_parameter["background_correction"] = background
  cell_detection_parameter["intensity_detection"]["measure"] = ["source"];
  cell_detection_parameter["shape_detection"]["threshold"] = cellshape #for 640 ch 500
  cell_detection_parameter["maxima_detection"]["shape"] = maximashape
  cell_detection_parameter["maxima_detection"]["save"] = False #DO NOT SAVE MAXIMA WTF

  processing_parameter = cells.default_cell_detection_processing_parameter.copy();
  processing_parameter.update(
        processes = 1, # 'serial', #multiple processes don't work on kepecs desktop bc of memory
        size_max = 10, #100, #35,
        size_min = 5,# 30, #30,
        optimization = False,
        optimization_fix = None,
        overlap  = 3, #32, #10,
        verbose = True
        )
  dstch2 = os.path.join(args.dst, "{0}_cells_{1}_raw.npy".format(args.brainname, args.channel2))
  print("\n\n STARTING CELL DETECTION FOR CHANNEL: {0}".format(args.channel2))       
  cells.detect_cells(ws.source("ch2"), dstch2,
                       cell_detection_parameter=cell_detection_parameter, 
                       processing_parameter=processing_parameter)
  
  #update
  # ws.update(cells_ch1 = dstch1, cells_ch2 = dstch2)
  #%% Filter cells
  
  thresholds = {
      "source" : None,
      "size"   : (20,None)
      }
  
  #channel 1
  cells.filter_cells(source = dstch1, 
                         sink = os.path.join(args.dst, "{0}_cells_{1}_filtered20.npy".format(args.brainname, args.channel1)), 
                         thresholds=thresholds);
  #channel 2
  cells.filter_cells(source = dstch2, 
                         sink = os.path.join(args.dst, "{0}_cells_{1}_filtered20.npy".format(args.brainname, args.channel2)), 
                         thresholds=thresholds);
  