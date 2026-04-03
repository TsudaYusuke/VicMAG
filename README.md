
<img width="200" src="https://github.com/TsudaYusuke/VicMAG/blob/main/regular.png">

# VicMAG   

VicMAG, an open-source tool for visualizing circular metagenome-assembled genomes highlighting bacterial virulence and antimicrobial resistance

## Prerequired
Python with the following libraries
- Pandas
- matplotlib
- Biopython
- Biotite
- pillow (PIL)

~~~
conda install pandas matplotlib -y
conda install -c conda-forge pillow biopython biotite -y
~~~

## Before run
 - extract cicular MAGs
 - annotate cMAGs by [DFAST](https://github.com/nigyta/dfast_core) (use option for ARG and VFG identification)
 - analyze cMAGs by geNomad, PlasFlow and CheckV (optional)


## Install


## Command

~~~ 
python vicmag.py --gbks /path_to_dir_containing_gbk_files --n_row 10  --outdir /path_to_output_dir
~~~

### Options
|command| |
-----|-----
|--gbks| directory containing multiple genbank files (required)|
|--n_row| the number of cMAGs in the first line (default:5)|
|--outdir| output directory|
|--force| remove exisitng output directory|
|||
|--plasflow| path to plasflow file|
|--checkv_qua| path to checkv quality file|
|--checkv_pro| path to checkv prophage file|
|--genomad_p| path to genomad summary_plasmid|
|--genomad_v| path to genomad summary_virus|
|||
|--plasmid_only| make map of plasmids only|
|--non_plasmid_only| make map of non_plasmids only|
|--v_a_only| make map of cMAGs containing vfgs or args|
|--virus_only| make map of cMAGs containing virus area|
|||
|--c_arg| color of antimicrobial resistance genes (default:red)|
|--c_vfg| color of virulence factor genes (default:green)|
|--c_cd| color of cds (default:lightgrey)|
|--c_vir| color of antimicrobial resistance genes (default:blue, alpha=0.3)|
|--c_non_p| color of antimicrobial resistance genes (default:azure)|




## Example
<img src="https://github.com/TsudaYusuke/VicMAG/blob/main/Example.png">

## Reference
Please cite the following article. 

VicMAG, an open-source tool for visualizing circular metagenome-assembled genomes highlighting bacterial virulence and antimicrobial resistance.
Yusuke Tsuda, Yasuhiro Tanizawa, Thi My Hanh Vu, Yosuke Nishimura, Masaki Shintani, Haruka Abe, Futoshi Hasebe, Ikuro Kasuga, Miki Nagao, Masato Suzuki.
bioRxiv 2026.03.31.714378; doi: https://doi.org/10.64898/2026.03.31.714378

Thank you.

