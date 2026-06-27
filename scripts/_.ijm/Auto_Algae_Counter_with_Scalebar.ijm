/*
 * SOP SCRIPT: UNIVERSAL ALGAE COUNTER (Scale Bar Eraser + .nd2 Support)
 * Lab Protocol: 2026 Update
 */

macro "SOP Universal Eraser [F2]" {
    // --- CONFIGURATION ---
    minSize = 50;          
    boxSize = 512; 
    threshMethod = "Triangle"; 
    
    // --- SETUP ---
    run("Close All");
    if (isOpen("ROI Manager")) { selectWindow("ROI Manager"); run("Close"); }
    print("\\Clear");
    
    inputDir = getDirectory("Select the INPUT Folder");
    outputDir = getDirectory("Select the OUTPUT Folder");
    
    list = getFileList(inputDir);
    Array.sort(list); 

    // Create CSV and write header
    csvPath = outputDir + "Algae_Results_Eraser.csv";
    File.saveString("Sample No,Conc (Cells/mL),Average,Counts...\n", csvPath);

    currentSample = "";
    counts = newArray(0); 
    imageCount = 0;

    for (i = 0; i < list.length; i++) {
        filename = list[i];
        ext = toLowerCase(filename);
        
        if (endsWith(ext, ".nd2") || endsWith(ext, ".tif") || endsWith(ext, ".tiff") || 
            endsWith(ext, ".jpg") || endsWith(ext, ".jpeg") || endsWith(ext, ".png")) {
            
            imageCount++;
            path = inputDir + filename;

            // 1. Parsing Sample ID
            dotIndex = lastIndexOf(filename, ".");
            nameNoExt = substring(filename, 0, dotIndex);
            if (indexOf(nameNoExt, "_") >= 0) {
                parts = split(nameNoExt, "_");
                sampleID = parts[0]; 
            } else {
                sampleID = nameNoExt;
            }

            // 2. Data Grouping
            if (sampleID != currentSample) {
                if (currentSample != "") { saveRowToCSV(currentSample, counts, csvPath); }
                currentSample = sampleID;
                counts = newArray(0); 
            }
            
            // 3. Open Image (Hybrid Logic)
            setBatchMode(false); // Must be OFF to see the image for erasing
            if (endsWith(ext, ".nd2")) {
                run("Bio-Formats Importer", "open=[" + path + "] color_mode=Default view=Hyperstack stack_order=XYCZT use_virtual_stack");
            } else {
                open(path);
            }
            
            originalTitle = getTitle();
            if (nSlices > 1) { run("Z Project...", "projection=[Max Intensity]"); originalTitle = getTitle(); }

            // 4. Auto-Crop
            makeRectangle(getWidth()/2 - boxSize/2, getHeight()/2 - boxSize/2, boxSize, boxSize);
            run("Crop");

            // 5. MANUAL INTERACTION: Scale Bar Eraser
            setTool("rectangle");
            // Beep to alert the user
            beep();
            waitForUser("Scale Bar Eraser", "Draw a box over the SCALE BAR or text.\nThen click OK.\n(If no scale bar, just click OK)");
            
            if (selectionType() != -1) {
                run("Set...", "value=0"); // Paint selection black
                run("Select None");
            }

            // 6. Analysis (Turn BatchMode back ON for speed)
            setBatchMode(true);
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
            
            // 7. Save Evidence
            selectWindow(originalTitle);
            run("8-bit");
            roiManager("Show All with labels");
            run("Labels...", "color=Cyan font=14 show use draw");
            run("Flatten"); 
            saveAs("Jpeg", outputDir + "Checked_" + nameNoExt + ".jpg");
            
            // 8. Cleanup
            run("Close All"); 
            setBatchMode(false); // Reset for next image eraser
        }
    }

    // Save final sample
    if (currentSample != "") { saveRowToCSV(currentSample, counts, csvPath); }
    
    if (imageCount == 0) {
        showMessage("Error", "No images found.");
    } else {
        showMessage("Success!", "Processed " + imageCount + " images.\nResults: " + csvPath);
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