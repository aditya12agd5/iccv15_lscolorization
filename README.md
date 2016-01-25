This code is implementation of the colorization algorithm published in proceedings of 
International Conference on Computer Vision (ICCV), 2015:

"Learning Large-Scale Automatic Image Colorization",
By: Aditya Deshpande, Jason Rock and David Forsyth

Copyright (c) 2016 University of Illinois at Urbana Champaign. 
All rights reserved.
 
Permission to use, copy, modify and distribute this software and its documentation for 
educational purpose is hereby granted without fee, provided that the above copyright 
notice and this permission notice appear in all copies of this software and that you 
do not sell the software.

THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND, EXPRESSED, IMPLIED 
OR OTHERWISE.


=== Acknowledgements ===

1) Includes Mean-shift code (modified by us to use Gaussian Kernel) by
*  Bryan Feldman 02/24/06
*  MeanShift first appears in
*  K. Funkunaga and L.D. Hosteler, "The Estimation of the Gradient of a
*  Density Function, with Applications in Pattern Recognition"

2) Includes fast EM implementation by
*  Patrick P. C. Tsui,
*  PAMI research group
*  Department of Electrical and Computer Engineering
*  University of Waterloo, 
*  March, 2006

*  Michael Boedigheimer
*  Amgen
*  Dept of Computational Biology
*  Thousand Oaks CA, 91320
*  Dec, 2005

=== Instructions === 

* The code was developed/tested for Matlab R2014b/2015a. 

* Run the wrapper script ColorizeMain as follows:

  ColorizeMain(train_dir, test_dir, output_dir);

  e.g. ColorizeMain('./data/beach_small/Train/', './data/beach_small/Test/', './data/ColorizeOut/');

  You can pass an optional last argument denoting the stage_num.

  stage_num = 1 => Train and Test
  stage_num = 2 => Only Test, this expects some *.mat files generated during training to be present in output_dir

* Sample data of beach scenes containing only 5 Train and 2 Test images is placed in the directory "data/beach_small/"

* You can download the large dataset separately from http://vision.cs.illinois.edu/projects/lscolor/

* If you use this code or data, cite the following paper:

  @inproceedings{LSColorization:ICCV15,
 	author    = {Aditya Deshpande and Jason Rock and David Forsyth},
	title     = {Learning Large-Scale Automatic Image Colorization},
	booktitle = {ICCV},
	year      = {2015}
  }

