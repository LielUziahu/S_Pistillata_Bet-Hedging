// GFP_Summary_Per_File_FIXED.ijm
// One summary row per ND2 file: File, Area, Mean, IntDen

inputDir = getDirectory("Choose folder with .nd2 files");
outputDir = inputDir;
outputCSV = outputDir + "GFP_thresholded_summary.csv";

// Write CSV header
File.saveString("File,Area,Mean,IntDen\n", outputCSV);

list = getFileList(inputDir);

// Set measurement options
run("Set Measurements...", "area mean integrated redirect=None decimal=3");

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".nd2")) {
        fullPath = inputDir + list[i];
        print("Processing: " + list[i]);

        // Open only the first series
        run("Bio-Formats Importer", "open=[" + fullPath + "] autoscale color_mode=Grayscale view=Hyperstack stack_order=XYCZT series=0");

        // If Z-stack, use middle slice
        slices = nSlices();
        if (slices > 1) {
            setSlice(round(slices/2));
        }

        // Convert to 8-bit for thresholding
        run("8-bit");

        // Apply Otsu threshold and convert to mask
        setAutoThreshold("Otsu");
        setOption("BlackBackground", false);
        run("Convert to Mask");

        // Analyze particles (only thresholded foreground)
        run("Analyze Particles...", "size=10-Infinity show=Nothing display summarize");

        // Read values from Summary table
        area = getResult("Area", 0);
        mean = getResult("Mean", 0);
        intDen = getResult("IntDen", 0);

        // Save results as CSV line
        File.append(list[i] + "," + area + "," + mean + "," + intDen + "\n", outputCSV);

        // Cleanup
        run("Close All");
        run("Clear Results");
    }
}

print("✅ Done! Results saved to:\n" + outputCSV);

