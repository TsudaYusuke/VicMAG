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
from logging import getLogger,INFO,DEBUG,StreamHandler,Formatter,FileHandler

parser = argparse.ArgumentParser('Options to run VicMAG')

parser.add_argument('--gbks',help='path to directory containing genbank files',required=True)
parser.add_argument('--outdir',help='output directory',default='./')
parser.add_argument('--n_row',help='number of cMAGs in the top row',default=10)

parser.add_argument('--plasflow',help='path to plasflow file',default='')
parser.add_argument('--checkv_qua',help='path to checkv quality file',default='')
parser.add_argument('--checkv_pro',help='path to checkv prophage file',default='')
parser.add_argument('--genomad_p',help='genomad summary_plasmid file',default='')
parser.add_argument('--genomad_v',help='genomad summary_virus file',default='')

parser.add_argument('--force',help='remove exisitng outdir',action='store_true')

parser.add_argument('--plasmid_only',help='make map of plasmids only',action='store_true')
parser.add_argument('--non_plasmid_only',help='make map of non_plasmids only',action='store_true')
parser.add_argument('--v_a_only',help='make map of cMAGs containing vfgs or args',action='store_true')
parser.add_argument('--virus_only',help='make map of cMAGs containing virus area',action='store_true')

parser.add_argument('--c_arg',help='color of antimicrobial resistance genes (default:red)',default='red')
parser.add_argument('--c_vfg',help='color of virulence factor genes (default:green)',default='green')
parser.add_argument('--c_cds',help='color of cds (default:lightgrey)',default='lightgrey')
parser.add_argument('--c_vir',help='color of virus genes (default:blue, alpha=0.3)',default='blue')
parser.add_argument('--c_non_p',help='color of not plasmid (default:azure)',default='azure')
parser.add_argument('--png_dpi', type=int, default=50, help='DPI for individual cMAG PNG files (default: 50)')
parser.add_argument('--tiff_dpi', type=int, default=300, help='DPI for final TIFF output (default: 300)')
parser.add_argument('--tiff', help='also save final cMAGS image as TIFF', action='store_true')

parser.add_argument('--max_vfg_labels',help='maximum number of VFG labels before collapsing nearby VFG labels (default:8)',type=int,default=8)
parser.add_argument('--vfg_group_angle',help='angular distance threshold in radians for grouping nearby VFG labels (default:0.15)',type=float,default=0.15)
parser.add_argument('--max_vfg_cluster_gap_bp', type=int, default=5000, help='maximum genomic distance in bp between neighboring VFGs to collapse into a cluster')

args = parser.parse_args()

if args.force:
    if os.path.isdir(args.outdir):
        shutil.rmtree(args.outdir)
    else:
        pass
else:
    pass

os.makedirs(args.outdir+'/tmp',exist_ok=True)

logger = getLogger(__name__)
logger.setLevel(INFO)

handler = StreamHandler()
handler.setLevel(INFO)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

file_handler = FileHandler(args.outdir+'/VicMAG.log' , 'w')
file_handler.setLevel(INFO)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

def custom_feature_formatter(feature,gene_color=biotite.colors["orange"],ARG_color=args.c_arg,VFG_color=args.c_vfg,CDS_color=args.c_cds):
    label = feature.qual.get("label")
    if feature.key == "CDS":
        return True, CDS_color, "black", None 
    elif feature.key == "rep_origin":
        return True, "blue", "black", None
    elif feature.key == "ARG":
        return True, ARG_color, 'black', None
    elif feature.key == 'hypothetical protein':
        return True, 'lightgrey', 'black', None
    elif feature.key == "VFG":
        return True, VFG_color, 'black', None
        
def sort_gb_by_length(gb):
    pd_gbs_len = pd.DataFrame([[i.name, len(i.seq)] for i in gb.values()],columns=['n','l']).sort_values('l')
    return pd_gbs_len['n'].to_list()

def add_margin_tate(img, sa1, height):
    haba = np.array(img).shape[1]
    result = Image.new(img.mode,(haba+round(haba*0.05)*2,height+sa1),(255,255,255))
    result.paste(img,(round(haba*0.05),round(sa1/2)))
    return result
    
def remove_margin(img):
    # Convert to RGB to avoid RGB/RGBA channel mismatch during concatenation.
    img = img.convert('RGB')
    bg = Image.new('RGB', img.size, (255, 255, 255))
    diff = ImageChops.difference(img, bg)
    croprange = diff.getbbox()

    # If the image is completely white or bbox cannot be detected, return the original image.
    if croprange is None:
        return img

    crop_img = img.crop(croprange)
    return crop_img
    
def return_checkv_data():
    extract_description = re.compile(r'(.*)\_(\d)\s(.*?)\/')
    checkv = []
    fasta = SeqIO.parse(args.outdir+'/tmp/fasta_checkv.fasta','fasta')
    for i in fasta:
        checkv.append(list(extract_description.search(i.description).groups()))
    
    checkv = pd.DataFrame(checkv)
    sum_checkv = pd.read_csv(args.outdir+'/tmp/sum_checkv.csv',index_col=0)
    for i in checkv[0].unique():
        sum_checkv.loc[i,'position'] = ','.join(checkv[checkv[0]==i][2].to_list())
    sum_checkv = sum_checkv[sum_checkv.index.str.endswith('c')]
    sum_checkv = sum_checkv[sum_checkv['provirus']=='Yes']
    
    sum_cir_comp_or_highq = sum_checkv[(sum_checkv['checkv_quality']=='High-quality') | (sum_checkv['checkv_quality']=='Complete')]
    
    return sum_cir_comp_or_highq



def draw_outer_arc(ax, x_start, x_end, r=1.001, linewidth=2.4):
    """
    Draw an outer black arc for a collapsed VFG cluster.
    The arc indicates the genomic region covered by the grouped VFGs.
    """
    if pd.isna(x_start) or pd.isna(x_end):
        return

    x_start = float(x_start)
    x_end = float(x_end)

    if not np.isfinite(x_start) or not np.isfinite(x_end):
        return

    two_pi = 2 * np.pi
    x_start = x_start % two_pi
    x_end = x_end % two_pi

    if x_end >= x_start:
        theta_arc = np.linspace(x_start, x_end, 150)
        r_arc = np.full_like(theta_arc, r)
        ax.plot(theta_arc, r_arc, c='black', linewidth=linewidth, solid_capstyle='round')
    else:
        theta_arc1 = np.linspace(x_start, two_pi, 100)
        theta_arc2 = np.linspace(0, x_end, 100)
        ax.plot(theta_arc1, np.full_like(theta_arc1, r), c='black', linewidth=linewidth, solid_capstyle='round')
        ax.plot(theta_arc2, np.full_like(theta_arc2, r), c='black', linewidth=linewidth, solid_capstyle='round')


def circular_midpoint(x_start, x_end):
    """
    Return the midpoint angle between x_start and x_end on a circle.
    """
    two_pi = 2 * np.pi
    x_start = float(x_start) % two_pi
    x_end = float(x_end) % two_pi

    if x_end >= x_start:
        return (x_start + x_end) / 2
    return ((x_start + x_end + two_pi) / 2) % two_pi




def format_vfg_cluster_label(text_series):
    """
    Format clustered VFG labels as:
        VFGs
        gene1
        gene2
        gene3

    Duplicated names are removed while preserving order.
    """
    genes = [str(x) for x in text_series]
    genes = list(dict.fromkeys(genes))
    return "<VFGs>\n" + "\n".join(genes)



def collapse_vfg_labels(label, max_vfg_labels=8, angle_gap=0.15, seq_len=None, max_vfg_cluster_gap_bp=5000):

    if len(label) == 0:
        return label

    label = label.sort_values('x').reset_index(drop=True)

    required_cols = {'kind', 'x_start', 'x_end', 'is_cluster'}
    if not required_cols.issubset(set(label.columns)):
        return label

    arg_label = label[label['kind'] == 'ARG'].copy()
    vfg_label = label[label['kind'] == 'VFG'].copy()

    # Collapse VFGs only for chromosome-scale sequences.
    # For plasmid/small contigs (<1,000,000 bp), keep all VFG labels individually.
    if seq_len is None or seq_len < 1_000_000:
        out = pd.concat([arg_label, vfg_label], ignore_index=True)
        return out.sort_values('x').reset_index(drop=True)

    if len(vfg_label) <= max_vfg_labels:
        out = pd.concat([arg_label, vfg_label], ignore_index=True)
        return out.sort_values('x').reset_index(drop=True)

    vfg_label = vfg_label.sort_values('x').reset_index(drop=True)

    groups = []
    current_group = [0]

    for n in range(1, len(vfg_label)):
        prev = vfg_label.loc[n - 1]
        curr = vfg_label.loc[n]

        if seq_len is not None:
            # Distance between the end of the previous VFG and the start of the current VFG.
            # Negative values mean overlapping features and should be grouped.
            prev_end_bp = float(prev['x_end']) / (2 * np.pi) * seq_len
            curr_start_bp = float(curr['x_start']) / (2 * np.pi) * seq_len
            gap_bp = curr_start_bp - prev_end_bp

            close_enough = gap_bp <= max_vfg_cluster_gap_bp
        else:
            gap_angle = curr['x'] - prev['x']
            close_enough = gap_angle < angle_gap

        if close_enough:
            current_group.append(n)
        else:
            groups.append(current_group)
            current_group = [n]

    groups.append(current_group)

    collapsed_vfg = []

    for g in groups:
        sub = vfg_label.loc[g]

        if len(sub) < 3:
            for _, row in sub.iterrows():
                collapsed_vfg.append([
                    row['x'],
                    1.02,
                    row['text'],
                    'VFG',
                    row['x_start'],
                    row['x_end'],
                    False
                ])
        else:
            group_text = format_vfg_cluster_label(sub['text'])

            cluster_start = sub['x_start'].min()
            cluster_end = sub['x_end'].max()
            cluster_mid = circular_midpoint(cluster_start, cluster_end)

            collapsed_vfg.append([
                cluster_mid,
                1.02,
                group_text,
                'VFG',
                cluster_start,
                cluster_end,
                True
            ])

    collapsed_vfg = pd.DataFrame(
        collapsed_vfg,
        columns=['x', 'y', 'text', 'kind', 'x_start', 'x_end', 'is_cluster']
    )

    out = pd.concat([arg_label, collapsed_vfg], ignore_index=True)
    out = out.sort_values('x').reset_index(drop=True)
    return out


def label_required_gap(prev_text, curr_text, rippo=1.0):
    """
    Estimate whether two neighboring labels need vertical staggering.

    The necessary angular gap depends on label length and visual map size.
    Larger circular maps have more physical space per radian, so the
    required angular gap is reduced using rippo.
    """
    prev_text = str(prev_text)
    curr_text = str(curr_text)

    max_len = max(len(prev_text), len(curr_text))
    min_len_text = min(len(prev_text), len(curr_text))

    # Short labels such as PmrF, filK, and YojI usually need less spacing.
    if min_len_text <= 4:
        required_gap = 0.055 + 0.004 * max_len
    else:
        required_gap = 0.070 + 0.006 * max_len

    # Scale by visual map size. For larger maps, the same angular distance
    # gives more actual separation, so reduce the threshold.
    visual_scale = 1 / np.sqrt(max(float(rippo), 0.1))
    visual_scale = min(max(visual_scale, 0.50), 1.20)

    required_gap *= visual_scale
    required_gap = min(max(required_gap, 0.055), 0.22)

    return required_gap


def label_y_step(prev_text, curr_text, target_len, rippo=1.0):
    """
    Vertical offset used only when labels are predicted to overlap.
    """
    prev_text = str(prev_text)
    curr_text = str(curr_text)

    if len(prev_text) <= 4 or len(curr_text) <= 4:
        step = 0.06
    else:
        step = 0.14 if target_len > 1_000_000 else 0.18

    # Slightly reduce step for larger maps.
    step *= min(max(1 / np.sqrt(max(float(rippo), 0.1)), 0.65), 1.0)

    return step



def label_in_strict_pi_zone(x):
    """
    Return True for angle regions where labels tend to overlap more easily:
      0.0-0.1 * pi
      0.9-1.1 * pi
      1.9-2.0 * pi
    """
    r = float(x) / np.pi
    return (0.0 <= r <= 0.1) or (0.9 <= r <= 1.1) or (1.9 <= r <= 2.0)


def zone_adjusted_gap_and_step(prev_x, curr_x, required_gap, y_step):
    """
    Make overlap avoidance stricter only in the specified pi-ranges.
    """
    if label_in_strict_pi_zone(prev_x) or label_in_strict_pi_zone(curr_x):
        required_gap *= 1.70
        y_step *= 1.30
    return required_gap, y_step



def nice_length_step(seq_len, target_n_ticks=5):
    """
    Choose a length tick interval in 1, 2, or 5 x 10^n units.
    Examples:
      20,000 bp    -> 5,000
      1,000,000 bp -> 200,000
    """
    seq_len = int(seq_len)
    if seq_len <= 0:
        return 1

    raw = seq_len / max(target_n_ticks, 1)
    exponent = int(np.floor(np.log10(raw)))
    base = 10 ** exponent
    fraction = raw / base

    if fraction <= 1:
        nice = 1
    elif fraction <= 2:
        nice = 2
    elif fraction <= 5:
        nice = 5
    else:
        nice = 10

    return int(nice * base)


def make_legend(arg,vfg,vir,npl):
    fig = plt.figure(figsize=(6,2 + 1.2 * sum([vir,npl])),layout='tight')
    ax = plt.axes()
    
    arrow_arg = patches.Polygon([(1,8.5),(2.5,8.5),(3,8),(2.5,7.5),(1,7.5)],closed=True,ec=None,fc=arg)
    arrow_vfg = patches.Polygon([(1,6.5),(2.5,6.5),(3,6),(2.5,5.5),(1,5.5)],closed=True,ec=None,fc=vfg)
    
    ax.add_patch(arrow_arg)
    ax.add_patch(arrow_vfg)
    
    ax.text(3.5,8,'Antimicrobial resistance gene',fontsize=20,va='center')
    ax.text(3.5,6,'Virurence factor gene',fontsize=20,va='center')
    
    if npl:
        if vir:
            plt.scatter(2,2,s=1500,marker='o',c='azure',edgecolor='black')
            ax.text(3.5,2,'Not plasmid',fontsize=20,va='center')
        else:
            plt.scatter(2,4,s=1500,marker='o',c='azure',edgecolor='black')
            ax.text(3.5,4,'Not plasmid',fontsize=20,va='center')

    if vir:    
        ax.plot([1,3],[4,4],c='black')
        ax.plot([1.5,2.5],[4,4],lw=10,c=args.c_vir,alpha=0.5)
        ax.text(3.5,4,'Virus area',fontsize=20,va='center')

    plt.axis('off')
    ax.set_xlim([0,20])
    ax.set_ylim([10-4.5*sum([vir,npl]),9])

    plt.savefig(args.outdir+'/legend.png')

def make_map(accession,show_name=True,show_length=True,store='image/'):
    pf = ''
    result_checkv = ''
    geno_p=''
    geno_v=''
    pd_gbs_l = pd.read_csv(args.outdir+'/summary_gbk.csv',index_col=0)
    
    if os.path.isfile(args.outdir+'/tmp/pf.csv'):
        pf = pd.read_csv(args.outdir+'/tmp/pf.csv')
        
    if os.path.isfile(args.outdir+'/tmp/geno_p.csv'):
        geno_p = pd.read_csv(args.outdir+'/tmp/geno_p.csv')
    
    if os.path.isfile(args.outdir+'/tmp/geno_v.csv'):
        geno_v = pd.read_csv(args.outdir+'/tmp/geno_v.csv')
        
    if os.path.isfile(args.outdir+'/tmp/sum_checkv.csv') and os.path.isfile(args.outdir+'/tmp/fasta_checkv.fasta'):
        result_checkv = return_checkv_data()
    
    target = accession
    s = target

    ARG = []
    VFG = []
    
    ann = []
    
    #for DFAST
    for i in target.features:
        if i.type == 'CDS':
            if i.location.strand > 0:
                if 'note' in i.qualifiers:
                    c_or_v = [s for s in i.qualifiers['note'] if 'similar to' in s]
                    if c_or_v:
                        if 'CARD' in c_or_v[0]:  
                            if 'gene' in i.qualifiers:
                                GENE = i.qualifiers['gene'][0]
                            else:
                                GENE = re.compile(r'similar to (.*) in CARD').search(c_or_v[0]).group(1)
                            CDS = seq.Feature(
                                 "ARG",
                                [seq.Location(i.location.start,i.location.end)],
                                {'product':GENE,'label':GENE}
                            )
                            ARG.append(GENE)
                            ann.append(CDS)
                            
                        elif 'VF_ID' in c_or_v[0]:                       # ARG > VFG if identified as both ARG and VFG
                            if 'gene' in i.qualifiers:
                                GENE = i.qualifiers['gene'][0]
                            else:
                                GENE = re.compile(r'similar to (.*) in VF_ID').search(c_or_v[0]).group(1)
                            CDS = seq.Feature(
                                 "VFG",
                                [seq.Location(i.location.start,i.location.end)],
                                {'product':GENE,'label':GENE}
                                )
                            VFG.append(GENE)
                            ann.append(CDS) 
                            
                        else:
                            CDS = seq.Feature(
                            "hypothetical protein",
                            [seq.Location(i.location.start,i.location.end)],
                            {'product':i.qualifiers['product'][0]}
                            )
                            ann.append(CDS)
                    else:
                        CDS = seq.Feature(
                        "CDS",
                        [seq.Location(i.location.start,i.location.end)],
                        {'product':i.qualifiers['product'][0],'label':i.qualifiers['product'][0]}
                    )
                    ann.append(CDS)
                else:        
                    CDS = seq.Feature(
                        "CDS",
                        [seq.Location(i.location.start,i.location.end)],
                        {'product':i.qualifiers['product'][0],'label':i.qualifiers['product'][0]}
                    )
                    ann.append(CDS)
            else:
                if 'note' in i.qualifiers:
                    c_or_v = [s for s in i.qualifiers['note'] if 'similar to' in s]
                    if c_or_v:
                        if 'CARD' in c_or_v[0]:
                            if 'gene' in i.qualifiers:
                                GENE = i.qualifiers['gene'][0]
                            else:
                                GENE = re.compile(r'similar to (.*) in CARD').search(c_or_v[0]).group(1)
                            CDS = seq.Feature(
                                 "ARG",
                                [seq.Location(i.location.start,i.location.end,seq.Location.Strand.REVERSE)],
                                {'product':GENE,'label':GENE}
                            )
                            ARG.append(GENE)
                            ann.append(CDS)
                            
                        elif 'VF_ID' in c_or_v[0]:
                            if 'gene' in i.qualifiers:
                                GENE = i.qualifiers['gene'][0]
                            else:
                                GENE = re.compile(r'similar to (.*) in VF_ID').search(c_or_v[0]).group(1)
                            CDS = seq.Feature(
                                 "VFG",
                                [seq.Location(i.location.start,i.location.end,seq.Location.Strand.REVERSE)],
                                {'product':GENE,'label':GENE}
                                )
                            VFG.append(GENE)
                            ann.append(CDS)
                            
                        else:
                            CDS = seq.Feature(
                            "hypothetical protein",
                            [seq.Location(i.location.start,i.location.end,seq.Location.Strand.REVERSE)],
                            {'product':i.qualifiers['product'][0]}
                            )
                            ann.append(CDS)
                    else:
                        CDS = seq.Feature(
                            "CDS",
                            [seq.Location(i.location.start,i.location.end,seq.Location.Strand.REVERSE)],
                            {'product':i.qualifiers['product'][0],'label':i.qualifiers['product'][0]}
                            )
                        ann.append(CDS)
                else:        
                    CDS = seq.Feature(
                        "CDS",
                        [seq.Location(i.location.start,i.location.end,seq.Location.Strand.REVERSE)],
                        {'product':i.qualifiers['product'][0],'label':i.qualifiers['product'][0]}
                        )
                    ann.append(CDS)
    if len(ARG)>0:
        pd_gbs_l.loc[target.name,'arg'] = ','.join(ARG)
    if len(VFG)>0:
        pd_gbs_l.loc[target.name,'vfg'] = ','.join(VFG)
    
    annotation =seq.Annotation(ann)
    
    rippo = (len(target.seq)**(1/3)) / (min_len**(1/3)) #adjust size
    
    fig = plt.figure(figsize=(7*rippo,7*rippo),tight_layout=True)
    ax = fig.add_subplot(111, projection="polar")
    
    #coloring chromosome from plasflow
    if type(pf) == type(pd.DataFrame()):
        if pf[pf['contig_name']==target.name]['label'].values[0].startswith('chromosome'):
            r = np.full(100,0.98)
            theta = np.linspace(0,2*np.pi,100)
            ax.fill(theta,r,args.c_non_p)
        elif pf[pf['contig_name']==target.name]['label'].values[0].startswith('plasmid'):
            pd_gbs_l.loc[target.name,'plasmid'] = 'yes'
        else:
            pass
            
    #coloring other from genomad
    if type(geno_p) == type(pd.DataFrame()):
        if target.name in list(geno_p['seq_name']):
            pd_gbs_l.loc[target.name,'plasmid'] = 'yes'
        else:
            r = np.full(100,0.98)
            theta = np.linspace(0,2*np.pi,100)
            ax.fill(theta,r,args.c_non_p)
    
    graphics.plot_plasmid_map(
        ax, annotation, plasmid_size=len(s), tick_step=len(s),
        feature_formatter=custom_feature_formatter,
        omit_oversized_labels=True,
        spacing=0
        )
        
    ax.set_axis_off()
    
    plt.ylim(0,1.60)
    
    label = pd.DataFrame(columns=['x','y','text','kind','x_start','x_end','is_cluster'])
    position_x = []
    position_y = []
    x = 0
    for i in annotation.get_features():
        loc = [loc for loc in i.locs][0]
        start = loc.first
        end = loc.last
        mean_loc = mean([start, end])

        x_pos = 2*np.pi*mean_loc/len(s)
        x_start = 2*np.pi*start/len(s)
        x_end = 2*np.pi*end/len(s)

        if i.key=='ARG':
            if i.qual['product'] != 'hypothetical protein':
                position_x.append(x_pos)
                position_y.append(0.95)
                label.loc[x,:] = [x_pos,1.1,i.qual['product'],'ARG',x_start,x_end,False]
                x += 1
            else:
                pass
        elif i.key=='VFG':
            if i.qual['product'] != 'hypothetical protein':
                position_x.append(x_pos)
                position_y.append(0.95)
                label.loc[x,:] = [x_pos,1.1,i.qual['product'],'VFG',x_start,x_end,False]
                x += 1
            else:
                pass
        else:
            pass
                
    # Collapse only VFG labels when VFG labels are too many.
    # ARG labels are kept individually.
    label = collapse_vfg_labels(
        label,
        max_vfg_labels=args.max_vfg_labels,
        angle_gap=args.vfg_group_angle,
        seq_len=len(s),
        max_vfg_cluster_gap_bp=args.max_vfg_cluster_gap_bp
    )
    label = label.sort_values('x').reset_index(drop=True)

    if len(label) == 0:
        position_y = []
    else:
        position_y = [1.02]

    for n in range(1, len(label)):
        prev_x = label.loc[n-1, 'x']
        curr_x = label.loc[n, 'x']
        gap = curr_x - prev_x

        prev_text = str(label.loc[n-1, 'text'])
        curr_text = str(label.loc[n, 'text'])

        # Keep the original base settings from label_fix.py,
        # but make the condition stricter only in specific pi-ranges.
        required_gap = label_required_gap(prev_text, curr_text, rippo=rippo)
        y_step = label_y_step(prev_text, curr_text, len(target), rippo=rippo)
        required_gap, y_step = zone_adjusted_gap_and_step(prev_x, curr_x, required_gap, y_step)

        if gap > required_gap:
            position_y.append(1.02)
        else:
            position_y.append(position_y[-1] + y_step)

    label['y'] = position_y


    if len(label) >= 2:
        first_x = label.loc[0, 'x']
        last_x = label.loc[len(label)-1, 'x']
        if (label_in_strict_pi_zone(first_x) and label_in_strict_pi_zone(last_x)):
            wrap_gap = (2 * np.pi - last_x) + first_x
            wrap_required_gap = label_required_gap(str(label.loc[len(label)-1, 'text']), str(label.loc[0, 'text']), rippo=rippo)
            wrap_y_step = label_y_step(str(label.loc[len(label)-1, 'text']), str(label.loc[0, 'text']), len(target), rippo=rippo)
            wrap_required_gap, wrap_y_step = zone_adjusted_gap_and_step(last_x, first_x, wrap_required_gap, wrap_y_step)

            if wrap_gap <= wrap_required_gap:
                label.loc[0, 'y'] = max(label.loc[0, 'y'], label.loc[len(label)-1, 'y'] + wrap_y_step)

    for i in label.index[::-1]:
        x0 = label.loc[i, 'x']
        y0 = label.loc[i, 'y']
        text = label.loc[i, 'text']
        is_cluster = bool(label.loc[i, 'is_cluster']) if 'is_cluster' in label.columns else False

        x_text = x0
        y_text = y0

        if x0 < np.pi / 2:
            ha = 'left'
            va = 'bottom'

        elif x0 < np.pi:
            ha = 'left'
            va = 'top'

        elif x0 < np.pi * 1.5:
            ha = 'right'
            va = 'top'

        else:
            ha = 'right'
            va = 'bottom'

        if is_cluster:
            x_start = label.loc[i, 'x_start']
            x_end = label.loc[i, 'x_end']

            # Keep the cluster label close to the outer arc.
            y_text = max(1.02, y0)

            # Draw the actual VFG cluster range as an outer black arc.
            draw_outer_arc(ax, x_start, x_end, r=1.001, linewidth=2.4)

            # Connect the label from just outside the arc.
            x_mid = circular_midpoint(x_start, x_end)
            ax.plot(
                [x_mid, x_text],
                [1.003, y_text],
                c='black',
                linewidth=0.8
            )
        else:
            ax.plot(
                [x0, x_text],
                [0.95, y_text],
                c='black',
                linewidth=0.8
            )

        ax.text(
            x_text,
            y_text,
            text,
            fontsize=7,
            bbox=dict(boxstyle="square,pad=0.15", ec='white', fc='white'),
            ha=ha,
            va=va
        )

    # show length
    show_len_haba = nice_length_step(len(s.seq), target_n_ticks=5)
    tick_positions = list(range(0, len(s.seq), show_len_haba))

    for i in tick_positions:
        x_num = 2*np.pi*i/len(s.seq)

        if x_num == 0:
            ax.plot([x_num,x_num],[1,0.80],c='black')
            ax.text(x_num,0.79,i,ha='center',va='top')
        elif x_num <np.pi*0.45:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='right',va='top')
        elif x_num <np.pi*0.55:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='right',va='center')
        elif x_num <np.pi*0.95:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='right',va='bottom')
        elif x_num < np.pi*1.05:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='center',va='bottom')
        elif x_num <np.pi*1.45:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='left',va='bottom')
        elif x_num <np.pi*1.55:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='left',va='center')
        elif x_num < np.pi*1.92:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='left',va='top')
        else:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.78,i,ha='right',va='top')

    #show virus area
    if type(result_checkv) == type(pd.DataFrame()):
        if target.name in result_checkv.index:
            prophage_area = result_checkv.loc[target.name,'position']
            for i in prophage_area.split(','):
                s_v, e_v = np.array(i.split('-')).astype(int)
                r_v = np.full(100,0.98)
                theta_v = np.linspace(2*np.pi*s_v/len(s),2*np.pi*e_v/len(s),100)
                ax.plot(theta_v,r_v,args.c_vir,linewidth=20,alpha=0.3)
            pd_gbs_l.loc[target.name,'virus'] = ','.join(prophage_area)
        else:
            pass
            
    if type(geno_v) == type(pd.DataFrame()):
        if list(geno_v[geno_v['seq_name'].str.contains(target.name)]['coordinates'].fillna(0)) != []:
            if list(geno_v[geno_v['seq_name'].str.contains(target.name)]['coordinates'].fillna(0)) != [0]:
                prophage_area = geno_v[geno_v['seq_name'].str.contains(target.name)]['coordinates'].values
                for i in prophage_area:
                    s_v, e_v = np.array(i.split('-')).astype(int)
                    r_v = np.full(100,0.98)
                    theta_v = np.linspace(2*np.pi*s_v/len(s),2*np.pi*e_v/len(s),100)
                    ax.plot(theta_v,r_v,'blue',linewidth=20,alpha=0.3)
                pd_gbs_l.loc[target.name,'virus'] = ','.join(prophage_area)
        else:
            pass
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.text(0,0.08,accession.id+'\n\n'+"{:,}".format(len(target))+' bp',fontsize=20*rippo/2,ha='center',va='top')
    
    plt.savefig(store+target.name+'.png',dpi=args.png_dpi)
    
    plt.close()
    
    pd_gbs_l.to_csv(args.outdir+'/summary_gbk.csv')

def page_main():
	pf = ''
	result_checkv = ''
	if os.path.isfile(args.outdir+'/tmp/pf.csv'):
		pf = pd.read_csv(args.outdir+'/tmp/pf.csv')
	
	if os.path.isfile(args.outdir+'/tmp/sum_checkv.csv') and os.path.isfile(args.outdir+'/tmp/fasta_checkv.fasta'):
		result_checkv = return_checkv_data()
	
	os.makedirs(args.outdir+'/tmp/image/',exist_ok=True)
	for i in gbs.keys():
		if not os.path.isfile(args.outdir+'/tmp/image/'+i+'.png'):
			logger.info('making each map: '+i)
			make_map(accession=gbs[i],store=args.outdir+'/tmp/image/')
		else:
			logger.info('There are already an image file. :'+i)
		
	exist_img_gbk = [Path(i).stem for i in os.listdir(args.outdir+'/tmp/image/')]
	not_exist = list(set(exist_img_gbk)-set(list(gbs.keys())))
	if not_exist:
		for i in not_exist:
			os.remove(args.outdir+'/tmp/image/'+i+'.png')
	else:
		pass
		
	logger.info('making maps')
	
	pd_gbs_l = pd.read_csv(args.outdir+'/summary_gbk.csv',index_col=0)
	
	if args.plasmid_only:
		pd_gbs_l = pd_gbs_l[pd_gbs_l['plasmid'] == 'yes']
	
	if args.non_plasmid_only:
		pd_gbs_l = pd_gbs_l[pd_gbs_l['plasmid'] != 'yes']
	
	if args.v_a_only:
		pd_gbs_l = pd_gbs_l[(~(pd_gbs_l['arg'].isnull())) | (~(pd_gbs_l['vfg'].isnull()))]
	
	if args.virus_only:
		pd_gbs_l = pd_gbs_l[~(pd_gbs_l['virus'].isnull())]
	
	imgs = []
	for i in pd_gbs_l.index[::-1]:
		imgs.append(remove_margin(Image.open(args.outdir+'/tmp/image/'+i+'.png')).convert('RGB'))
	
	# merge horizontal images in the top
	yoko = int(args.n_row)
	all_height = [np.array(i).shape[0] for i in imgs]
	height_1 = [np.array(i).shape[0] for i in imgs[:yoko]]
	img_1 = []
	for i in range(len(height_1)):
		sa_1 = max(height_1)-height_1[i]
		img_1.append(np.array(add_margin_tate(imgs[i],sa_1,height_1[i])))
	im_1 = np.concatenate(img_1,axis=1)
	max_width = im_1.shape[1]
	img_to_show = [im_1]
	
	# merge vertical images
	start = yoko
	x = 1
	ar = []
	all_width = [np.array(i).shape[1] for i in imgs]
	while start < len(all_width):
		while sum(all_width[start:start+x]) +2*round(sum(all_width[start:start+x]) *0.05) < max_width:
			if start+x == len(all_width)+1:
				ar.append(imgs[start:x+start-1])
				break
			x += 1
		else:
			ar.append(imgs[start:x+start-1])
		start = start + x -1 
		x = 1
	
	if len(ar) > 0:
		for k in ar:
			if len(k) == 0:
				continue
			height = [np.array(i.convert('RGB')).shape[0] for i in k]
			img_each = []
			for i in range(len(height)):
				sa_each = max(height)-height[i]
				img_each.append(np.array(add_margin_tate(k[i].convert('RGB'),sa_each,height[i])))
			ims = np.concatenate(img_each,axis=1)
			img_to_show.append(ims)
	else:
		# Keep img_to_show as a list.
		# If the number of cMAGs is <= n_row, there is only one horizontal row.
		img_to_show = [im_1]
	
	fig, ax = plt.subplots()
	ax.axis('off')
	
	# Pad each horizontal row to the same width before vertical concatenation.
	# Use the same dtype and channel number as the image array to avoid concat errors.
	img_yoko_width_each = [i.shape[1] for i in img_to_show]
	imgs_tate = []
	for i in img_to_show:
		sa_t = max(img_yoko_width_each)-i.shape[1]
		if sa_t > 0:
			pad = np.full((i.shape[0], sa_t, i.shape[2]), 255, dtype=i.dtype)
			i = np.concatenate([i, pad], axis=1)
		imgs_tate.append(i)
	im = np.concatenate(imgs_tate,axis=0)
	final_img = Image.fromarray(im.astype(np.uint8))
	final_img.save(args.outdir+'/cMAGS.png')
	if args.tiff:
		final_img.save(args.outdir+'/cMAGS.tiff', dpi=(args.tiff_dpi, args.tiff_dpi)) 
	
	logger.info('Done! See you!')

gbs_select = []
gbs = {}

if os.path.isdir(args.gbks):
	uploaded_files = os.listdir(args.gbks)
	if len(uploaded_files)>0:
		for uploaded_file in uploaded_files:
			if uploaded_file.endswith(('gb','gbk')):
				record = SeqIO.read(args.gbks+'/'+uploaded_file,'genbank')
				gbs[record.name] = record
			else:
				logger.warning('Unknown file:'+uploaded_file)
		pd_gbs_len = pd.DataFrame([[i.name, len(i.seq),None,None,None,None] for i in gbs.values()],columns=['n','l','plasmid','virus','arg','vfg']).sort_values('l')
		gbs_select = gbs_select + pd_gbs_len['n'].to_list()
		min_len = min(pd_gbs_len['l'])
		if os.path.isfile(args.outdir+'/summary_gbk.csv'):
			pass
		else:
			pd_gbs_len.to_csv(args.outdir+'/summary_gbk.csv',index=False)
	else:
		pass
		
f_plasmid = False
f_virus = False

if os.path.isfile(args.plasflow):
    pf = pd.read_table(args.plasflow,index_col=0)
    pf.to_csv(args.outdir+'/tmp/pf.csv')
    f_plasmid = True
else:
    logger.warning('No plasflow files')
    
if os.path.isfile(args.genomad_p):
    pf = pd.read_table(args.genomad_p,index_col=0)
    pf.to_csv(args.outdir+'/tmp/geno_p.csv')  
    f_plasmid = True  
else:
    logger.warning('No GenoVi plasmid file')
    
if os.path.isfile(args.genomad_v):
    pf = pd.read_table(args.genomad_v,index_col=0)
    pf.to_csv(args.outdir+'/tmp/geno_v.csv')  
    f_virus = True
else:
    logger.warning('No GenoVi virus file')

if os.path.isfile(args.checkv_qua) and os.path.isfile(args.checkv_pro):
            sum_checkv = pd.read_table(args.checkv_qua,index_col=0)
            sum_checkv.to_csv(args.outdir+'/tmp/sum_checkv.csv')
            
            fasta_checkv = SeqIO.parse(args.checkv_pro,'fasta')
            SeqIO.write(fasta_checkv,args.outdir+'/tmp/fasta_checkv.fasta','fasta')
            f_virus

make_legend(args.c_arg, args.c_vfg,f_virus,f_plasmid)

if len(gbs_select) > 1:
		page_main()
else:
	logger.warning('No file. Done')
