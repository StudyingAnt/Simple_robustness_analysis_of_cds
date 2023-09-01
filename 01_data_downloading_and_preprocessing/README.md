# Data downloading and preprocessing

## Downloading GENCODE 
To download GENCODE protein coding transcripts sequence, change `/path/to/download/` in the file `download_genecode.sh` to the path you want to download tha file and run the script. 

```console
./download_genecode.sh
```

## Preprocessing GENCODE
The file download with above script contains not only coding sequence (CDS) but also 5' untranslated region (5'UTR) and 3' untranslated region (3'UTR). 
In addition, several transcripts have deprecated coding sequences. This scripts extract CDS portion of valid transcripts and addtional 1 bp of 5' and 3' ends.

### Requirements
| Package | Version used in this study|
| --- | --- |
| Python | 3.9.16 |
| Biopython | 1.79 |

```console
python3 extract_cds_plus.py
```