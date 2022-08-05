


This directory stores generally used scripts and executables.


### ScriptManager-v0.13.jar
Downloaded by `/00_Download_and_preprocessing/1_download_files.sh`, this is the Java binary executable for ScriptManager that includes a collection of tools including TagPileup which is used to count tags and calculate coverage of samples around reference points.

### calculate_BED_ScoreRatio.pl
This script takes two input BED files (numerator BED file and denominator BED file) and writes an output BED entry for each gene in the files where the score column value is the ratio of values from the input files.
```
usage:		perl calculate_BED_ScoreRatio.pl  Numerator_BED_File	Denominator_BED_File	Output_BED
Example:	perl calculate_BED_ScoreRatio.pl numerator.bed denominator.bed ratio.bed
	BED information inherited from numerator BED file
```

### deduplicate_BED_coord_keep_highest_score.py
```
usage: deduplicate_BED_coord_keep_highest_score.py [-h] -i bedfile -o outfile

This script copies a BED coordinate file over to output while identifying duplicate coodinate entries and prints only the coordinate with the smallest absolute value in the score column to output.

optional arguments:
  -h, --help            show this help message and exit
  -i bedfile, --input bedfile
                        a BED file of coords to uniq (keep smallest absolute score)
  -o outfile, --output outfile
                        a BED file of deduplicated entries
```

### determine_closest_RefPoint_output_Both.pl

```
usage:		perl determine_closest_RefPoint_output_Both.pl	BED_File	Ref_BED_File	Output_File
Example:	perl determine_closest_RefPoint_output_Both.pl input.bed ref.bed output.bed
```

### filter_BED_by_list_ColumnSelect.pl
This script filters to keep or remove rows from any input BED (or tab-delimited) file based on the user-specified column values matching a user-specified list of strings.
```
usage:		perl filter_BED_by_list_ColumnSelect.pl	BED_File	List_File_Values	Column-Index (0-based)	keep/remove	Output_BED
Example:	perl filter_BED_by_list_ColumnSelect.pl input.bed ids.tab 4 keep output.bed
```

### filter_BED_by_string_ColumnSelect.pl
This script filters to keep or remove rows from any input BED (or tab-delimited) file based on the user-specified column values matching a user-specified string.
```
usage:		perl filter_BED_by_string_ColumnSelect.pl	BED_File	String value	Column-Index (0-based)	keep/remove	Output_BED
Example:	perl filter_BED_by_string_ColumnSelect.pl input.bed 01_RP 4 keep output.bed
```

### filter_BED_by_value_ColumnSelect.pl
This script filters to keep or remove rows from any input BED (or tab-delimited) file based on the user-specified column values being greater than or less than a numeric value.
```
usage:		perl filter_BED_by_value_ColumnSelect.pl	BED_File	Value (>=)	Column-Index (0-based)	keep/remove	Output_BED
Example:	perl filter_BED_by_value_ColumnSelect.pl input.bed 500 4 keep output.bed
```

### get_STM_from_Rossi_STable1.py
```
usage: get_STM_from_Rossi_STable1.py [-h] -i table -p peaks -o outbed [-s numbp]

This script processes STM occupancy reference out of Rossi STable1.

optional arguments:
  -h, --help            show this help message and exit
  -i table, --input table
                        The Rossi 2021 Supplemenatry Table 1 file as a tab-delimited flat text file
  -p peaks, --peaks peaks
                        Sua7 peaks as backup ref
  -o outbed, --output outbed
                        bed output
  -s numbp, --shift numbp
                        upstream bp shift when imputing default from TSS
```

### make_excel_composite_v2.py
This script takes a directory of composite (`*.out`) files and combines them into an excel spreadsheet with formatted metadata and convenient transformations.
```
usage: make_excel_composite_v2.py [-h] -i composite-dir -o outfile

This script takes a directory of composite (*.out) files and combines them into an excel spreadsheet.

optional arguments:
  -h, --help            show this help message and exit
  -i composite-dir, --input composite-dir
                        directory with all the composite data files
  -o outfile, --output outfile
                        output name to save workbook to
```

### shift_BED_center_v2.pl
This script shifts the coordinates of a BED file by a user-specified number of base pairs (+/- indicate upstream/downstream).
```
usage:		perl shift_BED_center.pl	BED_File	Shift(bp) [- = upstream, + = downstream]	Output_File
```

### sort_BED_by_Score_v2.pl
This script sorts the entries in a BED file by the score column (5th column, index-4) in either ascending or descending order.
```
usage:		perl sort_BED_by_Score_v2.pl	BED_File	desc/asc (desc = high->low, asc = low->high	Output_BED
Example:	perl sort_BED_by_Score_v2.pl input.bed asc output.bed
```

### sum_Col_CDT.pl
This script sums the columns of a CDT matrix file by column values (CDT to composite).
```
usage:		perl sum_Col_CDT.pl	Input_CDT_File	Output_TAB_File
Example:	perl sum_Col_CDT.pl input.cdt composite.out
```

### sum_Row_CDT.pl
This script sums the columns of a CDT matrix file by row values (CDT to occupancy).
```
usage:		perl sum_Row_CDT.pl	Input_CDT_File	Output_TAB_File
Example:	perl sum_Row_CDT.pl input.cdt ouput.tab
```

### update_BED_coord.pl
This script updates an input BED file with the coordinates of master BED file (lookup by identifier, 4th column, index-3).
```
usage:		perl update_BED_coord.pl	BED_File (to be updated)	BED_File (master list)	Output_BED
Example:	perl update_BED_coord.pl input.bed master.bed output.bed
```

### update_BED_score_with_TAB_score.pl
This script updates an input BED file's score column (5th column, index-4) with values from a 2-column tab-delimited input lookup ("id\tvalue").
```
usage:		perl update_BED_score_with_TAB_score.pl	BED_File (to be updated)	TAB_File (output from sum_Row_CDT.pl)	Output_BED
Example:	perl update_BED_score_with_TAB_score.pl input.bed ref.tab output.bed
```
