//different with v4, v5 will close all windows before processing the next image.
// Ask user to select input folder
imageFolder = getDirectory("Choose folder with .tif files: ");

// Create output folder for results
parentFolder = File.getParent(imageFolder);
folderName = File.getName(imageFolder);
resultsFolder = imageFolder + File.separator + folderName + "_Results";
File.makeDirectory(resultsFolder);

// Get list of all files in the image folder
fileList = getFileList(imageFolder);

// Loop through each file in the folder
for (i = 0; i < fileList.length; i++) {
    // Process only TIFF files (case-insensitive)
    if (endsWith(fileList[i].toLowerCase(), ".tif")) {
        // Define the current image and ROI paths
        imagePath = imageFolder + File.separator + fileList[i];
        roiPath = imageFolder + File.separator + replace(fileList[i], ".tif", "_rois.zip");
        resultsPath = resultsFolder + File.separator + replace(fileList[i], ".tif", ".csv");

        // Open the image
        open(imagePath);
        imageName = fileList[i];
        selectImage(imageName);

        // Load ROI file into ROI Manager (skip if ROI file doesn't exist)
        if (File.exists(roiPath)) {
            roiManager("Open", roiPath);

            // Select ROIs
            roiManager("Select All");
            roiManager("Show All");

            // Run BIOP JACoP analysis
            run("BIOP JACoP", 
                "channel_a=2 channel_b=3 "
                + "threshold_for_channel_a=MaxEntropy "
                + "threshold_for_channel_b=MaxEntropy "
                + "manual_threshold_a=0 manual_threshold_b=0 "
                + "get_pearsons get_spearmanrank get_manders get_overlap get_li_ica "
                + "costes_block_size=5 costes_number_of_shuffling=100"
            );

            // Save the results table
            saveAs("Results", resultsPath);

            // Manage ROIs
            roiManager("Select All");
            roiManager("Delete");

        } else {
            // Log missing ROI file
            print("ROI file not found for: " + imageName);
        }

        // Close the image after processing
        close();

        // Close all remaining windows (e.g., ROI Manager, results tables, etc.)
        run("Close All");
    }
}

// Print completion message
print("Batch processing completed! Results saved to: " + resultsFolder);
showMessage("Conversion done");

// Exit batch mode if it was enabled
setBatchMode(false);