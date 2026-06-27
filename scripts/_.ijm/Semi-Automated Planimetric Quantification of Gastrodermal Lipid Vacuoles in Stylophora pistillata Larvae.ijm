// =========================================================================
// 🎛️ CONTROL PANEL
// =========================================================================
var minDropletSize   = 5;
var blobSizeCutoff   = 300; // Increase this if lipids are misclassified as voids
var minCircularity   = 0.3;
var maxCircularity   = 1.00;
var blurRadius       = 1.5;
// =========================================================================
// 1. Setup: Directs the script to the folder with your raw OIR files and where to export CSVs.
macro "Batch Larvae Consumption Production v10.1" {
    inputDir = getDirectory("Select Folder Containing Larvae .oir Files");
    outputDir = getDirectory("Select Local Target Folder for Results");
    list = getFileList(inputDir);
    outputFile = outputDir + "Larvae_Consumption_Data_v10.1.csv";
    
// 2. CSV Initialization: Creates the CSV header file if it doesn't already exist.
    if (!File.exists(outputFile)) {
        File.saveString("File_Name,Total_Area,Endoderm_Area,Intact_Droplet_Area,Consumed_Blob_Area,Total_Vacuole_Area,Consumption_Index,Intact_Count,Blob_Count,Vacuole_Percentage_Total\n", outputFile);
    }
    
    setOption("BlackBackground", true);
// 3. Batch Loop: Cycles through every file in the directory.    
    for (i = 0; i < list.length; i++) {
    	// Only process files ending in .oir
        if (endsWith(list[i], ".oir") || endsWith(list[i], ".OIR")) {
            if (indexOf(File.openAsString(outputFile), list[i]) >= 0) continue; 
// 4. Import: Load the OIR file using the Bio-Formats library to preserve metadata.            
            run("Bio-Formats Importer", "open=[" + inputDir + list[i] + "] autoscale");
            mainTitle = getTitle(); 
            run("Split Channels");
            //Selecting Green channal for beter contrast
            selectWindow("C2-" + mainTitle); 
            rename("Working_Copy");
            close("\\Others");
// 5. Image Pre-processing: Convert to grayscale and apply Gaussian blur to prepare for segmentation.
            run("8-bit");
            

            setTool("polygon");
            waitForUser("Trace 1: Total Perimeter", "Trace outer perimeter.");
            run("Measure");
            totalArea = getResult("Area", 0);
            run("Clear Results");
            
            setTool("polygon");
            waitForUser("Trace 2: Endoderm", "Trace inner boundary.");
            run("Measure");
            endodermArea = getResult("Area", 0);
            run("Clear Results");
            
            // Masking background
            setBackgroundColor(0, 0, 0); 
            run("Clear Outside");
            
            roiManager("reset");
            roiManager("Add"); // Index 0: Endoderm
            run("Select None");
            
// 6. Thresholding (Manual adjustment gate)
            run("Duplicate...", "title=Base_Mask");
            run("Gaussian Blur...", "radius=" + blurRadius);
            
// 7.Force manual thresholding to solve "all red" background issues
            run("Threshold...");
            waitForUser("Verify Threshold", "Adjust slider so ONLY lipids/voids are red (and background is black), You may select diffrent threshold model (`Otsu` is recommended). DO NOT CLICK THE BUTTONS AT THE BOTTOM. Click OK to confirm.");
            run("Convert to Mask");
            
            run("Watershed");
            
// 8. Particle Analysis: Detects voids based on the size and circularity parameters set above.
            run("Set Measurements...", "area display");
            run("Analyze Particles...", "size=" + minDropletSize + "-Infinity circularity=" + minCircularity + "-" + maxCircularity + " display clear add");
            
            intactArea = 0; intactCount = 0;
            blobArea = 0; blobCount = 0;
            
// 9. Logic Loop: Categorizes results into 'Blobs' (irregular) or 'Droplets' (round).
            run("Clear Results");
            for (j = 0; j < roiManager("count"); j++) {
                roiManager("select", j);
                run("Measure");
                area = getResult("Area", nResults-1);
                
                if (area > blobSizeCutoff) {
                    blobArea += area;
                    blobCount++;
                    roiManager("Set Color", "magenta");
                    roiManager("Set Line Width", 3);
                    roiManager("Rename", "Void");
                } else {
                    intactArea += area;
                    intactCount++;
                    roiManager("Set Color", "cyan");
                    roiManager("Set Line Width", 1);
                    roiManager("Rename", "Drop");
                }
            }
            run("Clear Results");
            
// 10. Proof Generation: Saves a JPEG with ROI overlays to visually audit measurements.
            selectWindow("Working_Copy");
            roiManager("Show All without labels");
            run("Flatten");
            saveAs("Jpeg", outputDir + replace(list[i], ".oir", "_Consumption_Proof.jpg"));
            
// 11. Final Calculations: Computes the indices for export.
            totalVacuoleArea = intactArea + blobArea;
            
            consumptionIndex = 0;
            if (totalVacuoleArea > 0) {
                consumptionIndex = blobArea / totalVacuoleArea;
            }
            
            vacuolePercentTotal = 0;
            if (totalArea > 0) {
                vacuolePercentTotal = (totalVacuoleArea / totalArea) * 100;
            }
            
// 12. Logging: Appends all calculated metrics to the CSV file.
            File.append(list[i] + "," + totalArea + "," + endodermArea + "," + intactArea + "," + blobArea + "," + totalVacuoleArea + "," + consumptionIndex + "," + intactCount + "," + blobCount + "," + vacuolePercentTotal + "\n", outputFile);
            
// 13. Memory Cleanup: Resets ROIs and closes images to keep processing fast.
            roiManager("reset");
            run("Clear Results");
            run("Close All");
        }
    }
    showMessage("Success");
}