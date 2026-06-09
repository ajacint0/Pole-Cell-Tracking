# CAAX Embryo Tracking Protocol

## Use updated segmentation tool to segment nuclei up to timepoint where they cellularize
1. Pick a timepoint to determine when most of the pole cells cellularize
2. Obtain both ch 0 and ch 1 images of the movie you will be tracking up to the timepoint of cellularization
   - If you have a working CAAX StarDist model, evaluate on the movie and get the predictions to instead correct the segmentations
3. Put the ch 0 images in 'ch_0' folder in New_GUI, put ch 1 images in 'raw' folder in New_GUI
   - If you have model predictions, put those in 'seeds' folder in New_GUI
4. Run run_segmentation_tool.py after changing 'tp' variable to the timepoint you want to segment and 'ch' should be 1

## For timepoints after point of cellularization, use Cellpose for segmentation
1. Obtain ch 0 images of the movie you will be tracking after the timepoint of cellularization
2. 
