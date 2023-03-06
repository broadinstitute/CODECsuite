## Setup input for workflow
* sample_sheet.csv is used for `codex demux`, it requires three columns `SampleName,IndexBarcode1,IndexBarcode2`. Each row must have unique SampleName
* input.tsv store paths for fastq files and other pipeline inputs. The `sample` column in input.tsv has to match the `SampleName` column in sample_sheet.csv
