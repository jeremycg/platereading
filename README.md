platereading
============

r script for reading plates off a bioanalyser to QTL output.
Plates are fit to the buchanan curve by nls (there will be bad fits) and the parameters output for further analysis.

Takes input of a directory of directories of reads. The main directory must have one directory, labelled "plate xx", for each run. These directories must contain the .csv files given by the bioanalyser software, labelled xxx.csv, where xxx is the time in hours the reading was taken.
Additionally, the root directory must contain a "strainlist.csv" which has run, column, strain and temperature columns.

output is as .csv - the batchplates.R give "grouped.csv" which has all data grouped by line. This file is then displayed by the shiny app.
batchplates.Rmd gives the same output as "outputfits.csv", and "grouped.csv" finds means and sds based on a given residual cutoff.

The shiny app allows visualisation of data based on residual cutoffs - and provides a similar table as the batchplates.Rmd

todo:
Better integration of the batchplates.Rmd into shiny
More intelligent handling of bad fits
residuals based on number of time points
better handling of different data directories
