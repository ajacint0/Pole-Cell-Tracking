# CAAX Embryo Tracking Protocol
## Create virtual environment
1. Create environment with command ```python3 -m venv name_of_virtual_environment```
2. Activate the environment with ```source /mnt/home/YOUR_USERNAME/name_of_your_virtual_environment/bin/activate```
3. Install the python packages with ```pip install -r requirements.txt```
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
5. Save your segmentations in the format tp_{timepoint}_ch_1_seg.tif and transfer them to  /mnt/home/{user}/ceph/tracked_embryos/{movie_name}/gui_segs/. Create gui_segs if the folder doesn't exist

## For timepoints after point of cellularization, use Cellpose for segmentation
1. Obtain cropped 16/8bit ch 0 images of the movie you will be tracking and put in a folder in the movie name directory titled 'test_cellpose'
2. Make a new directory called 'cellpose' in the same directory as 'test_cellpose' is in
3. Edit cellpose_test.py and change the imwrite() function's path to point to 'cellpose' in your own directory
4. Also change varriable 'split_tp' to the timepoint of cellularization
5. Run cellpose_test.py

## Code for building tracking graph

1. Open match_nucleus_tp_new.py
   - Change 'user' variable to your ceph username
   - Change 'movie' variable to the name of the movie you will be tracking
   - Change 'split_tp' variable to the timepoint of cellularization
   - In the 'timevect' variable, set the range to the first and last timepoint you will be tracking
2. Run match_nucleus_tp_new.py
3. Open build_tree_from_match_files.py
   - Change 'user' variable to your ceph username
   - Change 'split_tp' variable to the timepoint of cellularization
   - Change 'movie' variable to the name of the movie you will be tracking
   - In the 'timevect' variable, set the range to the first and last timepoint you will be tracking
4. Run build_tree_from_match_files.py
  - This will create csv's which store the nodes and edges of the graph, edges denote connections between the same pole cells through time and nodes are each pole cell in a segmentation
5. Open visualize_edges.py
   - Change 'user' variable to your ceph username
   - Change 'movie' variable to the name of the movie you will be tracking
   - Set 'register' variable to True or False depending on if you want to view the early timepoint registered to the later one or not
   - Change 'split_tp' variable to the timepoint of cellularization
   - Changing the 'cutoff' variable to, for example, 20, will only highlight tracks that have at least 20 nodes connected through time
   - Changing 'limit_graph' variable to True will show all potential cells and tracks
   - Changing 'important_arr' variable to True will only show tracks that contain specific cells found at a certain timepoint, these cell labels are kept in a .npy file
   - Select timepoints to compare with the variables 'tp_early' and 'tp_late', they can be timepoints that are far apart such as 10 and 50

## Tracking Nuclei
1. Have tabs open for graph.csv, edges.csv, and the raw data + segmentations for the timepoints you are tracking on Fiji
   - Run visualize_edges.py taking note of the variables you set previously
   - An interactive window will pop up which will show centroids of pole cell nuclei differentiated by color.
        - Red corresponds to nuclei from the earlier timepoint
             - early nuclei which have lineages that trace back to nuclei in the 'important_arr' variable will be yellow
        - Blue corresponds to nuclei from the later timepoint
             - late nuclei which have lineages that trace back to nuclei in the 'important_arr' variable will be neon green
        - What you will want to do if you have a .npy file of important nuclei is to connect yellow nuclei with no links to blue nuclei, thus turning the blue nuclei green upon reloading the visualization
2. If you want to delete a connection between 2 nuclei, left click on the 2 nuclei that make the edge
3. If you want to add a connection, left click on the nuclei you want to connect
4. If you want to delete a nucleus, right click on it
5. To save the changes t
6. To add a nucleus to the segmentation, open edit_img.py and set the variable 'add_nucleus' to True enter and the zyx coordinates in the 'position' variable
7. To delete a nucleus from the segmentation, open edit_img.py and set the variable 'delete_nucleus' to True set the 'label' variable to the desired label
8. To split a nucleus in the segmentation, open edit_img.py and set the variable 'split_nucleus' to True and also set either 'cut_z', 'cut_y', or 'cut_x' to True. Edit either 'cutoffz', 'cutoffy', or 'cutoffx' to the frame you would like the cut to pass through
9. To combine nuclei in the segmentation, edit_img.py and set the variable combine_nuclei to True and set 'label' to the label the combined nuclei will have. Set 'false_label' to the label that will be merged into the other.
10. If splitting or adding a nucleus, you will have to add the new label into graph.csv in order for it to show up in the visualization, the same goes for deleting and merging nuclei but instead of adding, you will delete the label from graph.csv

<img width="1671" height="1670" alt="Screenshot from 2026-06-11 17-54-34" src="https://github.com/user-attachments/assets/56e1a32f-c30e-46bf-b4a5-4526420e219a" />
## Visualizing Tree
1. Open 'lineage_tree_new'
2. Change 'user' variable to your ceph username
3. Change 'movie' variable to the name of the movie you will be tracking
4. In the 'timevect' variable, set the range to the first and last timepoint you will be tracking
5. in the 'cuttoff' variable, enter the minimum number of nodes in a track to have it show up 
## Example Data
1. There is a folder with example segmentations and raw data in '/mnt/ceph/users/ajacinto/tracked_embryos/test/'.
2. Copy the 'test' folder into a directory called 'tracked_embryos' in your ceph
3. Start the protocol at 'Code for building tracking graph'
ADD HAYDEN ARRAY TO REPOSITORY
