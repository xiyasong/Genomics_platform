## Order of Running Scripts

To ensure proper execution and data processing flow, please run the scripts in the following order:

1. `Functions.R`
   - Description: A script contains several functions for faster processing two cohorts on the same way.

2. `Read_process_from_nodup4.R`
   - Description: Reads and processes data, utilizing functions defined in `Functions.R`.

3. `Whole_analysis_python_output.Rmd`
   - Description: A Rmd file to generate the figures in the manuscript.
  
4. `Database_structure_plot.R`
   - supplementary figure 1, gene panels relevant figures   
