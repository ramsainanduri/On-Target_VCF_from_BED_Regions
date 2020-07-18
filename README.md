# Extract On-Target and Off-Target VCF Records

Script to extract vcf records that falls in the given bed region. It also gives the vcf records which does not fall in the given bed regions (reverse or -r parameter).

## Requirements
* Python 3

## Inputs
The file inputs are 
* Original VCF file 
* BED file 

## Outputs
The file outputs are 
* On-Target VCF 

## Usage
```python
python On-Target_VCF_from_BED_Regions.py -i sample.vcf -b sample.bed -o output -r N
```
## Example
I have provided an example data set.
