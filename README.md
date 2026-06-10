# CAAX Embryo Tracking Protocol
## Prepare folder
1. Create directory 'tracked_embryos' in your ceph
2. In 'tracked_embryos', create another directory named after the movie you will be tracking

## Use updated segmentation tool to segment nuclei up to timepoint where they cellularize
1. Pick a timepoint to determine when most of the pole cells cellularize
2. Obtain both ch 0 and ch 1 cropped 16/8bit images of the movie you will be tracking up to the timepoint of cellularization
   - If you have a working CAAX StarDist model, evaluate on the movie and get the predictions to instead correct the segmentations
3. Put the ch 0 images in 'ch_0' folder in New_GUI, put ch 1 images in 'raw' folder in New_GUI
   - If you have model predictions, put those in 'seeds' folder in New_GUI
4. Run run_segmentation_tool.py after changing 'tp' variable to the timepoint you want to segment and 'ch' should be 1
5. Save your segmentations in the format tp_{timepoint}_ch_1_seg.tif

## For timepoints after point of cellularization, use Cellpose for segmentation
1. Obtain cropped 16/8bit ch 0 images of the movie you will be tracking after the timepoint of cellularization, put in a folder in the movie name directory titled 'test_cellpose'
2. Make a new directory called 'cellpose' in the same directory as 'test_cellpose' is in
3. Edit cellpose_test.py and change the imwrite() function's path to point to 'cellpose' in your own directory
4. Run cellpose_test.py

## Code for building tracking graph

1. Create directory called 'matches' in the same directory as 'cellpose' and 'test_cellpose' are in
2. Open match_nucleus_tp_new.py
   - Change occurances of 'ajacinto' to your ceph username
   - Change 'movie' variable to the name of the movie you will be tracking
   - Change 'split_tp' variable to the timepoint of cellularization
   - In the 'timevect' variable, set the range to the first and last timepoint you will be tracking
   - Change 'path_to_nuclear_segmentations' to the path where your saved segmentations from the tool are in
   - Change 'path_to_membrane_segmentations' to the path where your 'cellpose' folder is
   - Change 'path_for_saving_matches' to the path where your 'matches' folder is
4. Run match_nucleus_tp_new.py
5. Open build_tree_from_match_files.py
   - Change occurances of 'ajacinto' to your ceph username
   - Change 'split_tp' variable to the timepoint of cellularization
   - Change 'movie' variable to the name of the movie you will be tracking
   - In the 'timevect' variable, set the range to the first and last timepoint you will be tracking
