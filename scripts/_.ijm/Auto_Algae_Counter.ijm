/*
 * SOP SCRIPT: UNIVERSAL ALGAE COUNTER (ND2 + TIFF)
 * - Automatically detects .nd2 or .tif files
 * - Uses Bio-Formats for ND2, Standard Open for TIFF
 * - Crops to 1360x1360 (Middle Big Square)
 *change boxSize acording to the image pixel size:
 *if image is 1608x1680 set box to 1360
 *if image is 536x536 set box to 512
 */

macro "SOP Algae Universal [F4]" {
    // --- CONFIGURATION ---
    minSize = 50;           
    boxSize = 1360;  // Calculated size for Middle Big Square 
    threshMethod = "Triangle"; 
     
    // --- SETUP ---
    run("Close All");
    if (isOpen("ROI Manager")) { selectWindow("ROI Manager"); run("Close"); }
    print("\\Clear");
     
    inputDir = getDirectory("Select the INPUT Folder");
    outputDir = getDirectory("Select the OUTPUT Folder");
     
    list = getFileList(inputDir);
    Array.sort(list); 

    // Create CSV header
    csvPath = outputDir + "Algae_Results_Universal.csv";
    File.saveString("Sample No,Conc (Cells/mL),Average,Counts...\n", csvPath);

    currentSample = "";
    counts = newArray(0); 
    imageCount = 0;

    setBatchMode(true); 

    for (i = 0; i < list.length; i++) {
        filename = list[i];
        nameLower = toLowerCase(filename);
        
        // --- 1. CHECK FILE TYPE ---
        isND2 = endsWith(nameLower, ".nd2");
        isTIF = endsWith(nameLower, ".tif") || endsWith(nameLower, ".tiff");
        
        if (isND2 || isTIF) {
            imageCount++;
            showStatus("Processing: " + filename);

            // --- 2. PARSE SAMPLE ID ---
            // Remove extension intelligently
            nameNoExt = filename;
            if (isND2) { nameNoExt = replace(nameNoExt, ".nd2", ""); }
            if (isTIF) { 
                nameNoExt = replace(nameNoExt, ".tiff", ""); 
                nameNoExt = replace(nameNoExt, ".tif", ""); 
            }
            
            // Standardize Sample ID (split by underscore)
            if (indexOf(nameNoExt, "_") >= 0) {
                parts = split(nameNoExt, "_");
                sampleID = parts[0]; 
            } else {
                sampleID = nameNoExt;
            }

            // --- 3. DATA GROUPING ---
            if (sampleID != currentSample) {
                if (currentSample != "") { saveRowToCSV(currentSample, counts, csvPath); }
                currentSample = sampleID;
                counts = newArray(0); 
            }
            
            // --- 4. OPEN IMAGE (The Smart Part) ---
            path = inputDir + filename;
            
            if (isND2) {
                // Use Bio-Formats for ND2
                run("Bio-Formats Importer", "open=[" + path + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
            } else {
                // Use Standard Opener for TIFF
                open(path);
            }
            
            originalTitle = getTitle();
            
            // --- 5. CROP LOGIC (1360x1360 Center) ---
            makeRectangle(getWidth()/2 - boxSize/2, getHeight()/2 - boxSize/2, boxSize, boxSize);
            run("Crop");

            // --- 6. ANALYSIS ---
            run("Duplicate...", "title=Detection_Mask");
            run("8-bit");
            run("Subtract Background...", "rolling=50");
            setAutoThreshold(threshMethod + " dark");
            run("Convert to Mask");
            run("Watershed");
            
            // Count Particles
            roiManager("Reset");
            run("Analyze Particles...", "size=" + minSize + "-Infinity circularity=0.30-1.00 exclude add");
            
            thisCount = roiManager("count");
            counts = Array.concat(counts, thisCount);
            
            // --- 7. SAVE EVIDENCE ---
            selectWindow(originalTitle);
            run("8-bit"); 
            roiManager("Show All with labels");
            run("Labels...", "color=Cyan font=14 show use draw");
            run("Flatten"); 
            saveAs("Jpeg", outputDir + "Checked_" + nameNoExt + ".jpg");
            
            // Cleanup
            run("Close All"); 
        }
    }

    // Save final sample
    if (currentSample != "") { saveRowToCSV(currentSample, counts, csvPath); }
    
    setBatchMode(false);
    
    if (imageCount == 0) {
        showMessage("Error", "No .nd2 or .tif files found in the folder!");
    } else {
        showMessage("Success!", imageCount + " images processed.\nResults saved to: " + csvPath);
    }
}

function saveRowToCSV(sample, cArr, path) {
    sum = 0; n = cArr.length;
    for (k=0; k<n; k++) { sum = sum + cArr[k]; }
    avg = 0; if (n > 0) { avg = sum / n; }
    conc = avg * 10000; 
    line = sample + "," + conc + "," + avg;
    for (k=0; k<n; k++) { line = line + "," + cArr[k]; }
    File.append(line + "\n", path);
}