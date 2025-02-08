# VicMAG

VicMAG: visualizing circular metagenome-assembled genomes focused on bacterial virulence and antimicrobial resistance

>[!WARNING]
>This is alpha version. Now developing.

## Prerequired
Python with the following libraries
- Bioconda
- Biotite
- PIL

## Before run
 - extract cicular MAGs
 - annotate cMAGs by DFAST (use option for ARG and VFG identification)
 - analyze cMAGs by geNomad, PlasFlow and CheckV (optional)


## Install


## Command

~~~ 
python vicmag0201.py --dir  --yoko  --outdir
~~~

### Options
|command| |
-----|-----
|--dir| directory containing multiple genbank files (required)|
|--plasflow| plasflow file|
|--checkv_qua| checkv quality file|
|--checkv_pro| checkv prophage file|
|--genomad_p| genomad summary_plasmid|
|--genomad_v| genomad summary_virus|
|--n_row| the number of cMAGs in the first line (default:5)|
|--outdir| output directory|
|--force| remove exisitng output directory|

## Example

image

## Reference
Please 
