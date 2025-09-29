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
from PIL import Image
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
parser.add_argument('--c_vir',help='color of antimicrobial resistance genes (default:blue, alpha=0.3)',default='blue')
parser.add_argument('--c_non_p',help='color of antimicrobial resistance genes (default:azure)',default='azure')

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
    bg = Image.new('RGBA',img.size,(255, 255, 255))
    diff = ImageChops.difference(img,bg)
    croprange = diff.convert('RGB').getbbox()
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
    
    plt.ylim(0,1.5)
    
    label = pd.DataFrame(columns=['x','y','text'])
    position_x = []
    position_y = []
    x = 0
    for i in annotation.get_features():
        mean_loc = mean([[loc for loc in i.locs][0].first,[loc for loc in i.locs][0].last])
        if i.key=='ARG':
            if i.qual['product'] != 'hypothetical protein':
                position_x.append(2*np.pi*mean_loc/len(s))
                position_y.append(0.95)
                label.loc[x,:] = [2*np.pi*mean_loc/len(s),1.1,i.qual['product']]
                x += 1
            else:
                pass
        elif i.key=='VFG':
            if i.qual['product'] != 'hypothetical protein':
                position_x.append(2*np.pi*mean_loc/len(s))
                position_y.append(0.95)
                label.loc[x,:] = [2*np.pi*mean_loc/len(s),1.1,i.qual['product']]
                x += 1
            else:
                pass
        else:
            pass
                

    sa_between = [abs(j-i) for i, j in zip(np.sort(position_x)[:-1], np.sort(position_x)[1:])]
    position_y = [1.02]
    for i in sa_between:
        if i > 0.04:
            position_y.append(1.02)
        else:
            position_y.append(position_y[-1]+0.02)

    show_product = False
    if show_product:
        for i in sa_between:
            if i > 0.04:
                position_y.append(1.1)
            else:
                position_y.append(position_y[-1]+0.05)

    label = label.sort_values('x')
    
    label['y'] = position_y
    for i in label.index[::-1]:
        ax.plot([label.loc[i,'x'],label.loc[i,'x']],[0.95,label.loc[i,'y']],c='black')
        if label.loc[i,'x'] == 0:
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='center',
                   va='bottom')
        elif label.loc[i,'x'] <np.pi/2:
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='left',
                   va='bottom')
        elif label.loc[i,'x'] ==np.pi/2:
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='left',
                   va='center')
        elif label.loc[i,'x'] <np.pi:
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='left',
                   va='top')
        elif label.loc[i,'x'] == np.pi:
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='center',
                   va='bottom')
        elif label.loc[i,'x'] <np.pi*1.5:
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='right',
                   va='top')
        elif label.loc[i,'x'] ==np.pi*1.5:
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='right',
                   va='center')
        else:
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='right',
                   va='bottom')
                   
    # show length
    for i in range(0,len(s.seq),10**(len(str(len(s.seq)))-1)):
        x_num = 2*np.pi*i/len(s.seq)
        if x_num == 0:
            ax.plot([x_num,x_num],[1,0.80],c='black')
            ax.text(x_num,0.79,i,ha='center',va='top')
        elif x_num <np.pi/2:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='right',va='top')
        elif x_num ==np.pi/2:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='right',va='center')
        elif x_num <np.pi:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='right',va='bottom')
        elif x_num == np.pi:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='center',va='top')
        elif x_num <np.pi*1.5:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='left',va='bottom')
        elif x_num ==np.pi*1.5:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.79,i,ha='left',va='center')
        else:
            ax.plot([x_num,x_num],[0.82,0.80],c='black')
            ax.text(x_num,0.78,i,ha='left',va='top')
    
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
    
    plt.savefig(store+target.name+'.png',dpi=50)
    
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
		imgs.append(Image.open(args.outdir+'/tmp/image/'+i+'.png'))
	
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
			height = [np.array(i).shape[0]for i in k]
			img_each = []
			for i in range(len(height)):
			    sa_each = max(height)-height[i]
			    img_each.append(np.array(add_margin_tate(k[i],sa_each,height[i])))
			ims = np.concatenate(img_each,axis=1)
			img_to_show.append(ims)
	else:
		img_to_show = im_1
	
	fig, ax = plt.subplots()
	ax.axis('off')
	
	img_yoko_width_each = [i.shape[1] for i in img_to_show]
	imgs_tate = []
	for i in img_to_show:
		sa_t = max(img_yoko_width_each)-i.shape[1]
		imgs_tate.append(np.concatenate([i,np.full((i.shape[0],sa_t,i.shape[2]),255)],axis=1))
	im = np.concatenate(imgs_tate,axis=0)
	Image.fromarray(im.astype(np.uint8)).save(args.outdir+'/cMAGS.png') 
	
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

if os.path.isfile(args.checkv_qua) and os.path.file(args.checkv_pro):
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
