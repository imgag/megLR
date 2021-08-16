## Input Data

ONT Tools can analyse multiple samples in a single execution. All samples will be analysed with the same configuration. 

## Raw data folder

This pipeline is designed to run on a raw data folder of any ONT sequencing platform. The data can be available locally or, via a mounted network drive, still be on the machine. Before analysis, the reads are copied into the project folder. The raw data folders are not modified in any way. 

## File locations

ONT Tools normally requires a File of Filenames (FoF) to associate sample IDs and raw data locations. By default, the pipeline will look for a file named `sample_run_table.tsv` in the working directory (can be changed in options file). The FoF must be a Tab separated text file with either two or three columns:

1. Sample ID 
2. Data folder
3. Optional: Barcode. This column can be empty or omitted from the file. 

A single sample can have multiple data folders, this is used when multiple flowcells were used for the same sample, or if the sequencer created multiple runfolders because of interruptions. In this case there are multiple rows with the sampleID

SampleID | Folder  |   |
|  ---  |  ---  |  ---  |
| a0001    | /mnt/ontdata/run_20201205/      |       |
| a0001    | /mnt/ontdata/run_20201205_restarted/      |       |
| a0002    | /mnt/ontdata/run_20201206      |       |

The paths in the Folder column can be absolute or relative to the working directory. They will be searched recursively for all files ending with `.fastq`, `.fastq.gz`, `.fq` , `.fq.gz` and combined into a single compressed FastQ file in the sample directory. 

By default, subfolders with the name `failed` will be excluded from the search. This behaviour can be changed by setting `use_failed_reads = True` in the `config.yaml` file. 

### Running ONT Tools without FOFN

Of there is no valid FoF input, ONT Tools tries to gather Sample_IDs from existing folders. In this case, every Folder with a `Sample_` suffix will be searched for possible ONT analyses. 
