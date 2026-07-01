import os
import re
import shutil
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import biotite
import biotite.sequence as seq
import biotite.sequence.io.genbank as gb
import biotite.sequence.graphics as graphics
import biotite.database.entrez as entrez
from Bio import SeqIO
from PIL import Image, ImageChops, ImageOps
from pathlib import Path
from statistics import mean
from logging import getLogger, INFO, DEBUG, StreamHandler, Formatter, FileHandler


def main():
    parser = argparse.ArgumentParser('Options to run VicMAG')

    parser.add_argument('--gbks', help='path to directory containing genbank files', required=True)
    parser.add_argument('--outdir', help='output directory', default='./')
    parser.add_argument('--n_row', help='number of cMAGs in the top row', default=10)

    parser.add_argument('--plasflow', help='path to plasflow file', default='')
    parser.add_argument('--checkv_qua', help='path to checkv quality file', default='')
    parser.add_argument('--checkv_pro', help='path to checkv prophage file', default='')
    parser.add_argument('--genomad_p', help='genomad summary_plasmid file', default='')
    parser.add_argument('--genomad_v', help='genomad summary_virus file', default='')

    parser.add_argument('--force', help='remove exisitng outdir', action='store_true')

    parser.add_argument('--plasmid_only', help='make map of plasmids only', action='store_true')
    parser.add_argument('--non_plasmid_only', help='make map of non_plasmids only', action='store_true')
    parser.add_argument('--v_a_only', help='make map of cMAGs containing vfgs or args', action='store_true')
    parser.add_argument('--virus_only', help='make map of cMAGs containing virus area', action='store_true')

    parser.add_argument('--c_arg', help='color of antimicrobial resistance genes (default:red)', default='red')
    parser.add_argument('--c_vfg', help='color of virulence factor genes (default:green)', default='green')
    parser.add_argument('--c_cds', help='color of cds (default:lightgrey)', default='lightgrey')
    parser.add_argument('--c_vir', help='color of virus genes (default:blue, alpha=0.3)', default='blue')
    parser.add_argument('--c_non_p', help='color of not plasmid (default:azure)', default='azure')
    parser.add_argument('--png_dpi', type=int, default=50, help='DPI for individual cMAG PNG files (default: 50)')
    parser.add_argument('--tiff_dpi', type=int, default=300, help='DPI for final TIFF output (default: 300)')
    parser.add_argument('--tiff', help='also save final cMAGS image as TIFF', action='store_true')

    parser.add_argument('--max_vfg_labels', help='maximum number of VFG labels before collapsing nearby VFG labels (default:8)', type=int, default=8)
    parser.add_argument('--vfg_group_angle', help='angular distance threshold in radians for grouping nearby VFG labels (default:0.15)', type=float, default=0.15)
    parser.add_argument('--max_vfg_cluster_gap_bp', type=int, default=5000, help='maximum genomic distance in bp between neighboring VFGs to collapse into a cluster')

    args = parser.parse_args()

    # ここから下に、現在の line 54 以降の処理を入れる
    # if args.force:
    #     ...
    # make_legend(...)
    # if len(gbs_select) > 1:
    #     page_main()
    # else:
    #     logger.warning('No file. Done')


if __name__ == "__main__":
    main()