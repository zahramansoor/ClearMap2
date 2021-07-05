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
import sys, tifffile as tif
sys.path.append('/home/kepecs/python/ClearMap2/')

if __name__ == "__main__":
     
  #%############################################################################
  ### Initialization 
  ##############################################################################
  
  #% Initialize workspace
  
  from ClearMap.Environment import *  #analysis:ignore
  
  
  #directories and files
  directory = "/mnt/uncertainty/AA6-AK1a"    
    
  expression_raw      = "AA6-AK1a_561/original/AA6-AK1a_561-stitched_T001_Z<Z,I>_C01.tif"
  # "AA6-AK1a_640/AA6-AK1a_640-stitched_T001_Z<Z,I>_C01.tif"                        
  expression_auto     = "AA6-AK1a_488/AA6-AK1a_488-stitched_T001_Z<Z,I>_C01.tif"
  # directory = "/home/kepecs/Documents/2P_imaging/AA6-AK1a"    
    
  # expression_raw      = "cells_test/AA6-AK1a_640-stitched_T001_Z<Z,I>_C01.tif"            
  # expression_auto     = "561_cells_test/AA6-AK3b_561-stitched_T001_Z<Z,I>_C01.tif"            
  
  ws = wsp.Workspace("CellMap", directory=directory);
  ws.update(raw=expression_raw, autofluorescence=expression_auto)
  ws.debug = False
    
  resources_directory = settings.resources_path
    
  ws.info()
      
  resources_directory = settings.resources_path
  
  #%% Initialize alignment 
  
  #init atals and reference files
  annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(None)), orientation=(1,2,3),
      overwrite=False, verbose=True);
  
  #alignment parameter files    
  align_channels_affine_file   = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_affine_file  = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_bspline_file = io.join(resources_directory, 'Alignment/align_bspline.txt')
  
  
  #%%############################################################################
  ### Data conversion
  ############################################################################### 
  
  #%% Convet raw data to npy file     
               
  source = ws.source('raw');
  sink   = ws.filename('stitched')
  io.delete_file(sink)
  io.convert(source, sink, processes=None, verbose=True);
  
  
  #%%############################################################################
  ### Resampling and atlas alignment 
  ###############################################################################
        
  #%% Resample 
             
  resample_parameter = {
        "source_resolution" : (1.26, 1.26, 5.94), #z step might be 5 because less z planes
        "sink_resolution"   : (25,25,25),
        "orientation": (-3, -2, 1), #inverts old z (dorsal -> ventral) and flips x and z
        "processes" : None,
        "verbose" : True,             
        };
  #mod because of server priviledge issues
  sink = "/home/kepecs/Documents/AA6-AK1a_647_resampled.tif"   
  res.resample(ws.source('raw'), sink=sink, **resample_parameter) #sink=ws.filename('resampled')

  #%%
  p3d.plot(ws.filename('resampled'))
  
  #%% Resample autofluorescence
      
  resample_parameter_auto = {
        "source_resolution" : (1.26, 1.26, 5.94), #z step might be 5 because less z planes
        "sink_resolution"   : (25,25,25),
        "orientation": (-3, -2, 1), #inverts old z (dorsal -> ventral) and flips x and z
        "processes" : None,
        "verbose" : True,             
        };    
  #mod because of server priviledge issues
  sink = "/home/kepecs/Documents/AA6-AK1a_resampled.tif" 
  res.resample(ws.filename('autofluorescence'), sink=sink, **resample_parameter_auto) #ws.filename('resampled', postfix='autofluorescence')
  
  #p3d.plot([ws.filename('resampled'), ws.filename('resampled', postfix='autofluorescence')])
  
  #%% Aignment - resampled to autofluorescence
  
  # align the two channels
  align_channels_parameter = {            
      #moving and reference images
      "moving_image" : "/home/kepecs/Documents/AA6-AK1a_resampled.tif",
      "fixed_image"  : "/home/kepecs/Documents/AA6-AK1a_647_resampled.tif",
      
      #elastix parameter files for alignment
      "affine_parameter_file"  : align_channels_affine_file,
      "bspline_parameter_file" : None,
      
      #directory of the alig'/home/nicolas.renier/Documents/ClearMap_Ressources/Par0000affine.txt',nment result
      "result_directory" :  "/home/kepecs/Documents/AA6-AK1a_elastix_resampled_to_auto" 
      }; 
  
  elx.align(**align_channels_parameter);
  
  #%% Alignment - autoflourescence to reference
  
  # align autofluorescence to reference
  align_reference_parameter = {            
      #moving and reference images
      "moving_image" : "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_reference.tif", #whole brain
      "fixed_image"  : "/home/kepecs/Documents/AA6-AK1a_resampled.tif",
      
      #elastix parameter files for alignment
      "affine_parameter_file"  :  "/home/kepecs/python/ClearMap2/ClearMap/Resources/Alignment/align_affine.txt", 
      "bspline_parameter_file" :  "/home/kepecs/python/ClearMap2/ClearMap/Resources/Alignment/Order2_Par0000bspline.txt", #mods from brainpipe
      #directory of the alignment result
      "result_directory" :  "/home/kepecs/Documents/AA6-AK1a_elastix_auto_to_reference"
      };
  
  elx.align(**align_reference_parameter);
  
  
  #%%############################################################################
  ### Create test data
  ###############################################################################
  
  #%% Crop test data 
  
  #select sublice for testing the pipeline
  slicing = (slice(2000,2200),slice(2000,2200),slice(50,150));
  ws.create_debug('stitched', slicing=slicing);
  ws.debug = True; 
  
  #p3d.plot(ws.filename('stitched'))
    
  
  #%%############################################################################
  ### Cell detection
  ###############################################################################
  
  #%% Cell detection:
  from itertools import product
  #parameter sweep
  # bkshp = [(3,3)]
  # shpthres = [1900,2100,2300,2600]
  
  # # calculate number of iterations
  # tick = 0
  # for b, s in product(bkshp, shpthres):
  #    tick +=1
  # sys.stdout.write("\n\nNumber of iterations is {}:".format(tick))
  
  # for i in range(tick):
  #     bkshp_, shpthres_ = [xx for xx in product(bkshp, shpthres)][i] #parse out combinations
  #     print("\n")
  #     print("   iteration: {0}\n   background corr: {1}\n   shape detection: {2}".format(i,bkshp_[0], shpthres_))
  #     print("\n")
  cell_detection_parameter = cells.default_cell_detection_parameter.copy();
  cell_detection_parameter["illumination_correction"] = None;
  cell_detection_parameter["background_correction"] = None; #for 640 ch {"shape": (10,10), "form": "Disk"}; 
  cell_detection_parameter["intensity_detection"]["measure"] = ["source"];
  cell_detection_parameter["shape_detection"]["threshold"] = 2600 #for 640 ch 500
  # cell_detection_parameter["maxima_detection"]["shape"] = 10
  cell_detection_parameter["maxima_detection"]["save"] = False #DO NOT SAVE MAXIMA WTF
  # cell_detection_parameter["maxima_detection"]["save"] = ws.filename("cells", 
                                  # postfix="maxima_561_sweep_raw_bk{0}_shp{1}".format(bkshp_[0], shpthres_))
    # cell_detection_parameter['background_correction']['save'] = ws.filename("cells", 
    #                               postfix="background_561_sweep_raw_shp{0}".format(shpthres_))
    
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
  #old dest = ws.filename("cells", postfix="640_raw_bk10_shp500")
  dst = "/home/kepecs/Documents/AA6-AK1a_cells_561_raw_bkNone_shp2600.npy"    
  cells.detect_cells(ws.source("raw"), dst,
                       cell_detection_parameter=cell_detection_parameter, 
                       processing_parameter=processing_parameter)
  
  #%% visualization
  
  p3d.plot([[ws.filename('stitched'), ws.filename('cells', postfix='maxima')]])
  
  #%%
  coordinates = np.hstack([ws.source('cells', postfix='raw')[c][:,None] for c in 'xyz']);
  p = p3d.list_plot_3d(coordinates)
  p3d.plot_3d(ws.filename('stitched'), view=p, cmap=p3d.grays_alpha(alpha=1))
  
  
  #%% Cell statistics
  
  source = ws.source('cells', postfix='raw33_bk200')
  
  plt.figure(1); plt.clf();
  names = source.dtype.names;
  nx,ny = p3d.subplot_tiling(len(names));
  for i, name in enumerate(names):
    plt.subplot(nx, ny, i+1)
    plt.hist(source[name]);
    plt.title(name)
  plt.tight_layout();
  
  #%% Filter cells
  
  thresholds = {
      "source" : None,
      "size"   : (20,None)
      }
  
  #parameter sweep
  for i in range(tick):
      print(i)
      bkshp_, shpthres_ = [xx for xx in product(bkshp, shpthres)][i] #parse out combinations     
      cells.filter_cells(source = ws.filename("cells", postfix="561_sweep_raw_bk{0}_shp{1}".format(bkshp_[0], shpthres_)), 
                         sink = ws.filename("cells", postfix="561_sweep_raw_bk{0}_shp{1}_filtered_size20".format(bkshp_[0], shpthres_)), 
                         thresholds=thresholds);
  
  
  #%% Visualize
  
  coordinates = np.array([ws.source('cells', postfix='filtered')[c] for c in 'xyz']).T;
  p = p3d.list_plot_3d(coordinates, color=(1,0,0,0.5), size=10)
  p3d.plot_3d(ws.filename('stitched'), view=p, cmap=p3d.grays_alpha(alpha=1))
  
  
  #%%############################################################################
  ### Cell atlas alignment and annotation
  ###############################################################################
  
  #%% Cell alignment
  
  source = io.as_source("/home/kepecs/Documents/AA6-AK1a_cells_640_raw_bk10_shp500.npy")
  imgfl = "/mnt/uncertainty/AA6-AK1a/AA6-AK1a_561/original/"
  imgpth = [os.path.join(imgfl, xx) for xx in os.listdir(imgfl)]
  y,x = tif.imread(imgpth[0]).shape
  z = len(imgpth)
  rs = "/home/kepecs/Documents/AA6-AK1a_resampled.tif" 
  def transformation(coordinates):
    coordinates = res.resample_points(
                    coordinates, sink=None, orientation=(-3, -2, 1), 
                    source_shape=(x,y,z), 
                    sink_shape=io.shape(rs));
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory="/home/kepecs/Documents/AA6-AK1a_elastix_resampled_to_auto", 
                    binary=True, indices=False);
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory="/home/kepecs/Documents/AA6-AK1a_elastix_auto_to_reference",
                    binary=True, indices=False);
        
    return coordinates;
    
  
  coordinates = np.array([source[c] for c in 'xyz']).T;
  
  coordinates_transformed = transformation(coordinates);
  
  # Cell annotation
  
  label = ano.label_points(coordinates_transformed, key='order');
  names = ano.convert_label(label, key='order', value='name');
  
  # Save results
  
  coordinates_transformed.dtype=[(t,float) for t in ('xt','yt','zt')]
  label = np.array(label, dtype=[('order', int)]);
  names = np.array(names, dtype=[('name' , 'U256')])
  
  import numpy.lib.recfunctions as rfn
  cells_data = rfn.merge_arrays([source[:], coordinates_transformed, label, names], flatten=True, usemask=False)
  
  io.write("/home/kepecs/Documents/AA6-AK1a_cells.npy", cells_data)
  
  #############################################################################
  ### Cell csv generation for external analysis
  ###############################################################################
  
  # CSV export
  
  source = np.load("/home/kepecs/Documents/AA6-AK1a_cells.npy")
  header = ', '.join([h[0] for h in source.dtype.names]);
  np.savetxt("/home/kepecs/Documents/AA6-AK1a_cells.csv", source[:], header=header, delimiter=',', fmt='%s')
  
  #%% ClearMap 1.0 export
  
  source = ws.source('cells');
  
  clearmap1_format = {'points' : ['x', 'y', 'z'], 
                      'points_transformed' : ['xt', 'yt', 'zt'],
                      'intensities' : ['source', 'dog', 'background', 'size']}
  
  for filename, names in clearmap1_format.items():
    sink = ws.filename('cells', postfix=['ClearMap1', filename]);
    data = np.array([source[name] if name in source.dtype.names else np.full(source.shape[0], np.nan) for name in names]);
    io.write(sink, data);
  
  
  #%%############################################################################
  ### Voxelization - cell density
  ###############################################################################
  
  source = np.load("/home/kepecs/Documents/AA6-AK1a_points/posttransformed_zyx_voxels.npy")
  
  coordinates = np.array([[xx[2],xx[1],xx[0]] for xx in source]);#xyz
  intensities = source['source'];
  
  # Unweighted 
  
  #whole brain
  annotation_file = "/home/kepecs/python/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_annotation.tif"
  voxelization_parameter = dict(
        shape = io.shape(annotation_file), 
        dtype = None, 
        weights = None,
        method = 'sphere', 
        radius = (7,7,7), 
        kernel = None, 
        processes = None, 
        verbose = True
      )
  
  vox.voxelize(coordinates, sink="/home/kepecs/Documents/AA6-AK1a_density_counts.tif", **voxelization_parameter);
  
  
  #%% 
  
  p3d.plot(ws.filename('density', postfix='counts'))
  
  
  #%% Weighted 
  
  voxelization_parameter = dict(
        shape = io.shape(annotation_file),
        dtype = None, 
        weights = intensities,
        method = 'sphere', 
        radius = (7,7,7), 
        kernel = None, 
        processes = None, 
        verbose = True
      )
  
  vox.voxelize(coordinates, sink=ws.filename('density', postfix='intensities'), **voxelization_parameter);
  
  #%%
  
  p3d.plot(ws.filename('density', postfix='intensities'))
