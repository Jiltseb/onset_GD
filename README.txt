Copyright (C) 2018 Speech and Music Technology Lab, Indian Institute of Technology Madras
                
This file is part of GD based onset detection(onset_GD).                  
onset_GD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or   
(at your option) any later version.
                
This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.                             

You should have received a copy of the GNU General Public License. If not, see <http://www.gnu.org/licenses/>

This folder contains the main matlab program onset_GD.m which calls the binary "WordSegmentWithSilenceRemoval" to perform GD-based smoothing, the parameters of which can be adjusted in the configuration file: "fe-words.base_ref". 
The audio file should be placed inside the folder "test_here". The placement of the file is essential. 
The binary WordSegmentWithSielnceRemoval and the script extrema.m must be in the same folder.
This is the single resolution program, hence wsF is a parameter. 
The MLF output file (HTK format) is stored in this folder. Use mirex standards for the comparison with the ground truth.

More details about the work is discussed at a later work in https://sites.google.com/site/percussiononsetdetection/

If you use this code, please cite our paper at NCC 2015: https://ieeexplore.ieee.org/abstract/document/7084897/

P. A. Manoj Kumar, J. Sebastian and H. A. Murthy, "Musical onset detection on carnatic percussion instruments," 2015 Twenty First National Conference on Communications (NCC), Mumbai, 2015, pp. 1-6, doi: 10.1109/NCC.2015.7084897.



