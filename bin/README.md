


This directory stores generally used scripts and executables.


### ScriptManager-v0.13.jar
Downloaded by `/00_Download_and_preprocessing/1_download_files.sh`, this is the Java binary executable for ScriptManager that includes a collection of tools including TagPileup which is used to count tags and calculate coverage of samples around reference points.

### calculate_BED_ScoreRatio.pl
This script takes two input BED files (numerator BED file and denominator BED file) and writes an output BED entry for each gene in the files where the score column value is the ratio of values from the input files.

### deduplicate_BED_coord_keep_highest_score.py
### determine_closest_RefPoint_output_Both.pl
### filter_BED_by_list_ColumnSelect.pl
This script filters to keep or remove rows from any input BED (or tab-delimited) file based on the user-specified column values matching a user-specified list of strings.

### filter_BED_by_string_ColumnSelect.pl
This script filters to keep or remove rows from any input BED (or tab-delimited) file based on the user-specified column values matching a user-specified string.

### filter_BED_by_value_ColumnSelect.pl
This script filters to keep or remove rows from any input BED (or tab-delimited) file based on the user-specified column values being greater than or less than a numeric value.

### get_STM_from_Rossi_STable1.py
### shift_BED_center_v2.pl
This script shifts the coordinates of a BED file by a user-specified number of base pairs (+/- indicate upstream/downstream).

### sort_BED_by_Score_v2.pl
This script sorts the entries in a BED file by the score column (5th column, index-4) in either ascending or descending order.

### sum_Col_CDT.pl
This script sums the columns of a CDT matrix file by column values (CDT to composite).

### sum_Row_CDT.pl
This script sums the columns of a CDT matrix file by row values (CDT to occupancy).

### update_BED_coord.pl
This script updates an input BED file with the coordinates of master BED file (lookup by identifier, 4th column, index-3).

### update_BED_score_with_TAB_score.pl
This script updates an input BED file's score column (5th column, index-4) with values from a 2-column tab-delimited input lookup ("id\tvalue").
