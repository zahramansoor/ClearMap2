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
from scipy.io import loadmat
import sys, tifffile as tif, os, itertools, shutil, numpy as np
import argparse
sys.path.append('/home/kepecs/python/ClearMap2/')

def transformix_command_line_call(src, dst, transformfile):
    '''Wrapper Function to call transformix using the commandline, this can be time consuming
    
    Inputs
    -------------------
    src = volume path for transformation
    dst = folder to save file
    transformfile = final transform file from elastix registration
    
    '''
    from subprocess import check_output
    print ('Running transformix, this can take some time....\n')
    #sp.call(['transformix', '-in', src, '-out', dst, '-tp', transformfile])
    call = 'transformix -in {} -out {} -tp {}'.format(src, dst, transformfile)
    print(check_output(call, shell=True))
    print('Past transformix command line Call')      
            
    return

def change_transform_parameter_initial_transform(fl, initialtrnsformpth):
    '''
    (InitialTransformParametersFileName "NoInitialTransform")
    initialtrnsformpth = 'NoInitialTransform' or 'pth/to/transform.0.txt'
    '''
    fl1 = fl[:-5]+'0000.txt'
    with open(fl, 'r') as f, open(fl1, 'w') as new:
            for line in f:
                new.write('(InitialTransformParametersFileName "{}")\n'.format(initialtrnsformpth)) if 'InitialTransformParametersFileName' in line else new.write(line)
    #overwrite original transform file
    shutil.move(fl1, fl)
    return

def modify_transform_files(transformfiles, dst):
    """Function to copy over transform files, modify paths in case registration was done on the cluster, and tether them together
    
        Inputs
    ---------
    transformfiles = 
        list of all elastix transform files used, and in order of the original transform****
    """
    #new
    ntransformfiles = [os.path.join(dst, "order{}_{}".format(i,os.path.basename(xx))) for i,xx in enumerate(transformfiles)]    
    #copy files over
    [shutil.copy(xx, ntransformfiles[i]) for i,xx in enumerate(transformfiles)]   
    #modify each with the path
    for i,pth in enumerate(ntransformfiles):
        #skip first
        if i!=0:
            #read
            with open(pth, "r") as fl:
                lines = fl.readlines()
                fl.close()
            #copy
            nlines = lines
            #iterate
            for ii, line in enumerate(lines):
                if "(InitialTransformParametersFileName" in line:
                    nlines[ii] = "(InitialTransformParametersFileName {})\n".format(ntransformfiles[i-1])
            #write
            with open(pth, "w") as fl:
                for nline in lines:
                    fl.write(str(nline))
                fl.close()
    return ntransformfiles

def point_transformix(pretransform_text_file, transformfile, dst):
    """apply elastix transform to points      
    Inputs
    -------------
    pretransform_text_file = list of points that already have resizing transform
    transformfile = elastix transform file
    dst = folder
    
    Returns
    ---------------
    trnsfrm_out_file = pth to file containing post transformix points
    
    """
    sys.stdout.write("\n***********Starting Transformix***********")
    from subprocess import check_output
    #set paths    
    trnsfrm_out_file = os.path.join(dst, "outputpoints.txt")
    #run transformix point transform
    call = "transformix -def {} -out {} -tp {}".format(pretransform_text_file, dst, transformfile)
    print(check_output(call, shell=True))
    sys.stdout.write("\n   Transformix File Generated: {}".format(trnsfrm_out_file)); sys.stdout.flush()
    return trnsfrm_out_file

def create_text_file_for_elastix(src, dst):
    """
    Inputs
    ---------
    src = numpy file consiting of nx3 (ZYX points)
    dst = folder location to write points
    """
    print("This function assumes ZYX centers...")
    #setup folder
    if not os.path.exists(dst): os.mkdir(dst)
    #create txt file, with elastix header, then populate points
    pth=os.path.join(dst, "zyx_points_pretransform.txt")
    #load
    if type(src) == np.ndarray:
        arr = src
    else:
        arr = np.load(src) if src[-3:] == "npy" else loadmat(src)["cell_centers_orig_coord"]
    #convert
    stringtowrite = "\n".join(["\n".join(["{} {} {}".format(i[2], i[1], i[0])]) for i in arr]) ####this step converts from zyx to xyz*****
    #write file
    sys.stdout.write("writing centers to transfomix input points text file..."); sys.stdout.flush()
    with open(pth, "w+") as fl:
        fl.write("index\n{}\n".format(len(arr)))    
        fl.write(stringtowrite)
        fl.close()
    sys.stdout.write("...done writing centers\n"); sys.stdout.flush()
    return pth

def unpack_pnts(points_file, dst):
    """
    function to take elastix point transform file and return anatomical locations of those points
    
    Here elastix uses the xyz convention rather than the zyx numpy convention
    
    Inputs
    -----------
    points_file = post_transformed file, XYZ
    
    Returns
    -----------
    dst_fl = path to numpy array, ZYX
    
    """   
    #####inputs 
    assert type(points_file)==str
    point_or_index = 'OutputPoint'
    #get points
    with open(points_file, "r") as f:                
        lines=f.readlines()
        f.close()
    #####populate post-transformed array of contour centers
    sys.stdout.write("\n\n{} points detected\n\n".format(len(lines)))
    arr=np.empty((len(lines), 3))    
    for i in range(len(lines)):        
        arr[i,...]=lines[i].split()[lines[i].split().index(point_or_index)+3:lines[i].split().index(point_or_index)+6] #x,y,z            
    #optional save out of points
    dst_fl = os.path.join(dst, "posttransformed_zyx_voxels.npy")
    np.save(dst_fl, np.asarray([(z,y,x) for x,y,z in arr]))    
    #check to see if any points where found
    print("output array shape {}".format(arr.shape))
        
    return dst_fl
def parsecmdline():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "brainname", help="Brain name")
    parser.add_argument(
        "zstep", help="Z step size (in um)",
        type=float)
    parser.add_argument(
        "--autochannel", help="Autofluorescence channel", type=int, default=488)
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
        type=tuple, default=(3, -2, 1))
    parser.add_argument(
        "--dst", help="Save destination",
        type=str, default="/home/kepecs/Documents/")
    return parser.parse_args()

def fillargs(args):
    """kepecs lab experiment specific params
    for Shujing's cfos dataset, July 2021"""
    args.src = "/mnt/uncertainty"
    args.background640 = (10,10)
    args.cellshape640 = 500
    args.maximashape640 = 5
    args.background561 = None
    args.cellshape561 = 2600
    args.maximashape561 = 10
    if not os.path.exists(args.dst): os.mkdir(args.dst)
    
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
  str_auto = os.listdir(os.path.join(directory, "{0}_{1}".format(args.brainname, args.autochannel)))[98][:-12]
  str_ch1 = os.listdir(os.path.join(directory, "{0}_{1}".format(args.brainname, args.channel1)))[89][:-12]
  expression_ch1      = "{0}_{1}/{2}Z<Z,I>_C01.tif".format(args.brainname, args.channel1, str_ch1)
  expression_auto     = "{0}_{1}/{2}Z<Z,I>_C01.tif".format(args.brainname, args.autochannel, str_auto)

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
  align_channels_affine_file   = io.join(resources_directory, "Alignment/align_affine.txt")
  align_reference_affine_file  = io.join(resources_directory, "Alignment/align_affine.txt")
  align_reference_bspline_file = io.join(resources_directory, "Alignment/Order2_Par0000bspline.txt")
  
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
  res.resample(ws.source("ch1"), sink=sink_ch1, **resample_parameter) 
  if isinstance(args.channel2,int):
      sink_ch2 = os.path.join(args.dst, args.brainname+"_{}_resampled.tif".format(args.channel2))   
      res.resample(ws.source("ch2"), sink=sink_ch2, **resample_parameter) 

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
  res.resample(ws.filename("autofluorescence"), sink=sink_auto, **resample_parameter_auto) 
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
  
  elx.align(**align_channels_parameter);
  
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
      
      elx.align(**align_channels_parameter);
  
  
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
  
  elx.align(**align_reference_parameter);
  

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

#%% Point transform - zd edits
  #channel 1 & 2
  if isinstance(args.channel1,int) and isinstance(args.channel2, int):
      channels = [(1,args.channel1), (2,args.channel2)]
  elif isinstance(args.channel1,int) and not isinstance(args.channel2, int):
      channels = [(1,args.channel1)]
  for chnum, channel in channels:    
       cells = os.path.join(args.dst, "{0}_cells_{1}_filtered20.npy".format(args.brainname, channel))
       pnts = np.load(cells)
       pnts = np.array([pnts[c] for c in 'zyx']).T
       #for resize dimensions
       downsized = tif.imread(sink_ch1) #sagittal
       zd,yd,xd = downsized.shape #sagittal
       #reorient pnts
       pnts_sag = np.array([[xx[2],xx[1],xx[0]] for xx in pnts]) #from xyz
       #get full size dims
       stitched = os.path.join(directory, "{0}_{1}".format(args.brainname, channel))
       y,z = tif.imread(os.path.join(stitched, os.listdir(stitched)[56])).shape #sagittal
       x = len([xx for xx in os.listdir(stitched) if ".tif" in xx]) #sagittal
       f = ((zd/z),(yd/y),(xd/x))
       downsized_pnts_sag = np.array([[xx[0]*f[0],xx[1]*f[1],xx[2]*f[2]] for xx in pnts_sag]).astype(int)
       #map cells first
       cell=np.zeros((zd,yd,xd)) 
       for pnt in downsized_pnts_sag :
           z,y,x=pnt
           cell[z,y,x] = 1
       #flip x and y
       cell_oriented = np.flip(cell, axis=0)    
       #get coordinates again
       downsized_pnts_sag_oriented = np.nonzero(cell_oriented)
       downsized_pnts_sag_oriented = np.array([downsized_pnts_sag_oriented[0],downsized_pnts_sag_oriented[1],downsized_pnts_sag_oriented[2]]).T
       #transform
       #make into transformix-friendly text file
       transformed_dst = os.path.join(args.dst, "{0}_{1}_points".format(args.brainname, channel))
       if not os.path.exists(transformed_dst): os.mkdir(transformed_dst)
       pretransform_text_file = create_text_file_for_elastix(downsized_pnts_sag_oriented, transformed_dst)
       #get transform files
       transformfiles = []
       lstch1 = [os.path.join(os.path.join(args.dst, args.brainname+"_elastix_resampled_to_auto_ch{0}".format(chnum)), xx)
                 for xx in os.listdir(os.path.join(args.dst, args.brainname+"_elastix_resampled_to_auto_ch{0}".format(chnum))) if "TransformParameters" in xx]; lstch1.sort()
       lstauto = [os.path.join(os.path.join(args.dst, args.brainname+"_elastix_auto_to_reference" ), xx) 
                  for xx in os.listdir(os.path.join(args.dst, args.brainname+"_elastix_auto_to_reference" )) if "TransformParameters" in xx]; lstauto.sort()
       transformfiles.append(lstch1)
       transformfiles.append(lstauto)
       transformfiles = list(itertools.chain.from_iterable(transformfiles))
      
       #copy over elastix files
       transformfiles = modify_transform_files(transformfiles, transformed_dst) 
       change_transform_parameter_initial_transform(transformfiles[0], 'NoInitialTransform')
       #run transformix on points
       points_file = point_transformix(pretransform_text_file, transformfiles[-1], transformed_dst)
       #convert registered points into structure counts
       converted_points = unpack_pnts(points_file, transformed_dst)
 
