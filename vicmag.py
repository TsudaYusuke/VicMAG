from Bio import SeqIO
import os
from io import StringIO
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.io.genbank as gb
import biotite.sequence.graphics as graphics
import biotite.database.entrez as entrez
import biotite
import re
from PIL import Image
import subprocess
import math
from statistics import mean
import shutil
import tempfile
import argparse

def custom_feature_formatter(feature,gene_color=biotite.colors["orange"],AMG_color='red',VFG_color='green'):
    # AddGene stores the feature label in the '\label' qualifier
    label = feature.qual.get("label")
    if feature.key == "CDS":
        return True, 'lightgrey', "black", None   #gene_colorをlightgreyに変更
    elif feature.key == "rep_origin":
        return True, "blue", "black", None
    elif feature.key == "AMG":
        return True, AMG_color, 'black', None
    elif feature.key == 'hypothetical protein':
        return True, 'lightgrey', 'black', None
    elif feature.key == "VFG":
        return True, VFG_color, 'black', None
        
def remove_hypo(annotation):
    new_ann = []
    for i in annotation:
        if i.key != 'hypothetical protein':
            new_ann.append(i)
    return seq.Annotation(new_ann)
    
def return_checkv_data():
    extract_description = re.compile(r'(.*)\_(\d)\s(.*?)\/')
    checkv = []
    fasta = SeqIO.parse('tmp/fasta_checkv.fasta','fasta')
    for i in fasta:
        checkv.append(list(extract_description.search(i.description).groups()))
    
    checkv = pd.DataFrame(checkv)
    sum_checkv = pd.read_csv('tmp/sum_checkv.csv',index_col=0)
    for i in checkv[0].unique():
        sum_checkv.loc[i,'position'] = ','.join(checkv[checkv[0]==i][2].to_list())
    sum_checkv = sum_checkv[sum_checkv.index.str.endswith('c')]
    sum_checkv = sum_checkv[sum_checkv['provirus']=='Yes']
    
    sum_cir_comp_or_highq = sum_checkv[(sum_checkv['checkv_quality']=='High-quality') | (sum_checkv['checkv_quality']=='Complete')]
    
    return sum_cir_comp_or_highq

def make_map(accession,show_name=True,show_length=True,store='image/'):
    pf = ''
    result_checkv = ''
    geno_p=''
    geno_v=''
    
    if os.path.isfile('tmp/pf.csv'):
        pf = pd.read_csv('tmp/pf.csv')
        
    if os.path.isfile('tmp/geno_p.csv'):
        geno_p = pd.read_csv('tmp/geno_p.csv')
    
    if os.path.isfile('tmp/geno_v.csv'):
        geno_v = pd.read_csv('tmp/geno_v.csv')
        
    if os.path.isfile('tmp/sum_checkv.csv') and os.path.isfile('tmp/fasta_checkv.fasta'):
        result_checkv = return_checkv_data()
    
    #geNomad対応
    #dir内のファイル参照plasmidとvirus
    #csv読み取って、表示する
    
    target = accession
    s = target

    AMG = []
    
    ann = []
    
    
    #DFAST対応、VFG対応
    for i in target.features:
        if i.type == 'CDS':
            if i.location.strand > 0:
                if 'note' in i.qualifiers:
                    for j in i.qualifiers['note']:
                        if 'ARG' in j:              #dfastはCARDと書かれるので変更する
                            if 'gene' in i.qualifiers:
                                GENE = i.qualifiers['gene'][0]
                            else:
                                GENE = ''  #i.qualifiers['product'][0]
                            CDS = seq.Feature(
                                 "AMG",
                                [seq.Location(i.location.start,i.location.end)],
                                {'product':GENE,'label':GENE}
                            )
                            AMG.append(GENE)
                            ann.append(CDS)
                            #break
                        elif 'VFDB' in j:                       # CARDとVFDBが両方存在する時CARD優先
                            if 'gene' in i.qualifiers:
                                GENE = i.qualifiers['gene'][0]
                            else:
                                GENE = ''#i.qualifiers['product'][0]
                            CDS = seq.Feature(
                                 "VFG",
                                [seq.Location(i.location.start,i.location.end)],
                                {'product':GENE,'label':GENE}
                                )
                            AMG.append(GENE)
                            ann.append(CDS) #理想的にはidentityで載せるかどうか判定を入れるべき
                            #break
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
                if 'note' in i.qualifiers:
                    for j in i.qualifiers['note']:
                        if 'CARD' in j:
                            if 'gene' in i.qualifiers:
                                GENE = i.qualifiers['gene'][0]
                            else:
                                GENE = ''#i.qualifiers['product'][0]
                            CDS = seq.Feature(
                                 "AMG",
                                [seq.Location(i.location.start,i.location.end,seq.Location.Strand.REVERSE)],
                                {'product':GENE,'label':GENE}
                            )
                            AMG.append(GENE)
                            ann.append(CDS)
                            #break
                        elif 'VFDB' in j:
                            if 'gene' in i.qualifiers:
                                GENE = i.qualifiers['gene'][0]
                            else:
                                GENE = ''#i.qualifiers['product'][0]
                            CDS = seq.Feature(
                                 "VFG",
                                [seq.Location(i.location.start,i.location.end,seq.Location.Strand.REVERSE)],
                                {'product':GENE,'label':GENE}
                                )
                            AMG.append(GENE)
                            ann.append(CDS)
                            #break
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
                    
    annotation =seq.Annotation(ann)
    
    rippo = math.pow(len(target.seq),1/3) / math.pow(25893,1/3) #立方比とる
    
    fig = plt.figure(figsize=(7*rippo,7*rippo),tight_layout=True)
    ax = fig.add_subplot(111, projection="polar")
    
    #coloring chromosome
    if type(pf) == type(pd.DataFrame()):
        if pf[pf['contig_name']==target.name]['label'].values[0].startswith('chromosome'):
            r = np.full(100,0.98)
            theta = np.linspace(0,2*np.pi,100)
            ax.fill(theta,r,'azure')
        elif pf[pf['contig_name']==target.name]['label'].values[0].startswith('plasmid'):
            pass
        else:
            pass
    
    #coloring other from genomad
    if type(geno_p) == type(pd.DataFrame()):
        if geno_p.loc[geno_p['seq_name']==target.name]['plasmid_score'].values<0.95:
            r = np.full(100,0.98)
            theta = np.linspace(0,2*np.pi,100)
            ax.fill(theta,r,'azure')
        else:
            pass
    
    graphics.plot_plasmid_map(
        ax, annotation, plasmid_size=len(s), tick_step=len(s),
        feature_formatter=custom_feature_formatter,
        omit_oversized_labels=True,
        spacing=0
        )
        
    ax.set_axis_off()
    
    plt.ylim(0,1.5)
    label = pd.DataFrame(columns=['x','y','text'])
    texts = []
    position_x = []
    position_y = []
    x = 0
    for i in annotation.get_features():
        mean_loc = mean([[loc for loc in i.locs][0].first,[loc for loc in i.locs][0].last])
        if i.key=='AMG':
            if i.qual['product'] != 'hypothetical protein':
                position_x.append(2*np.pi*mean_loc/len(s))
                position_y.append(0.95)
                label.loc[x,:] = [2*np.pi*mean_loc/len(s),1.1,i.qual['product']]
                x += 1
    
    ticks = ax.get_xticks()
    labels = ax.get_xticklabels()
    sa_between = [abs(j-i) for i, j in zip(np.sort(position_x)[:-1], np.sort(position_x)[1:])]
    position_y = [1.1]
    for i in sa_between:
        if i > 0.04:
            position_y.append(1.1)
        else:
            position_y.append(position_y[-1]+0.05)

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
        if label.loc[i,'x'] <np.pi:
            ax.plot([label.loc[i,'x'],label.loc[i,'x']],[0.95,label.loc[i,'y']],c='black')
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='left')
        else:
            ax.plot([label.loc[i,'x'],label.loc[i,'x']],[0.95,label.loc[i,'y']],c='black')
            ax.text(label.loc[i,'x'],label.loc[i,'y'],label.loc[i,'text'],bbox=dict(boxstyle="square",
                   ec='white',
                   fc='white',
                   ),
                   ha='right')

    for i in range(0,len(s.seq),10**(len(str(len(s.seq)))-1)):
        x_num = 2*np.pi*i/len(s.seq)
        if x_num == 0:
            ax.text(x_num,1,i,ha='center',va='bottom')
        elif x_num <np.pi/2:
            ax.plot([x_num,x_num],[1,1.01],c='black')
            ax.text(x_num,1.02,i,ha='left',va='bottom')
        elif x_num ==np.pi/2:
            ax.plot([x_num,x_num],[1,1.01],c='black')
            ax.text(x_num,1.02,i,ha='left',va='center')
        elif x_num <np.pi:
            ax.plot([x_num,x_num],[1,1.01],c='black')
            ax.text(x_num,1.02,i,ha='left',va='top')
        elif x_num == np.pi:
            ax.plot([x_num,x_num],[1,1.01],c='black')
            ax.text(x_num,1.02,i,ha='center',va='top')
        elif x_num <np.pi*1.5:
            ax.plot([x_num,x_num],[1,1.01],c='black')
            ax.text(x_num,1.02,i,ha='right',va='top')
        elif x_num ==np.pi*1.5:
            ax.plot([x_num,x_num],[1,1.01],c='black')
            ax.text(x_num,1.02,i,ha='right',va='center')
        else:
            ax.plot([x_num,x_num],[1,1.01],c='black')
            ax.text(x_num,1.02,i,ha='right',va='bottom')
    
    #show virus area
    if type(result_checkv) == type(pd.DataFrame()):
        if target.name in result_checkv.index:
            prophage_area = result_checkv.loc[target.name,'position']
            for i in prophage_area.split(','):
                s_v, e_v = np.array(i.split('-')).astype(int)
                r_v = np.full(100,0.98)
                theta_v = np.linspace(2*np.pi*s_v/len(s),2*np.pi*e_v/len(s),100)
                ax.plot(theta_v,r_v,'blue',linewidth=20,alpha=0.3)
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
        else:
            pass
    
    ax.text(0,0,accession.id,fontsize=20*rippo/2,ha='center',va='top')#20*rippo +'\n\n'+"{:,}".format(len(target))+' bp'
    plt.savefig(store+target.name+'.png',dpi=100)
    
def sort_gb_by_length(gb):
    pd_gbs_len = pd.DataFrame([[i.name, len(i.seq)] for i in gb.values()],columns=['n','l']).sort_values('l')
    return pd_gbs_len['n'].to_list()

def add_margin(img, sa1, width):
    result = Image.new(img.mode,(width+sa1,width+sa1),(255,255,255))
    result.paste(img,(round(sa1/2),round(sa1/2)))
    return result

def page_main():
	pf = ''
	result_checkv = ''
	if os.path.isfile('tmp/pf.csv'):
		pf = pd.read_csv('tmp/pf.csv')
	
	if os.path.isfile('tmp/sum_checkv.csv') and os.path.isfile('tmp/fasta_checkv.fasta'):
		result_checkv = return_checkv_data()
		

	os.makedirs('tmp/image/',exist_ok=True)
	for i in gbs.keys():
		if not os.path.isfile('tmp/image/'+i+'.png'):
			make_map(accession=gbs[i],store='tmp/image/')
		else:
			pass
		
	exist_img_gbk = [Path(i).stem for i in os.listdir('tmp/image/')]
	not_exist = list(set(exist_img_gbk)-set(list(gbs.keys())))
	if not_exist:
		for i in not_exist:
			os.remove('tmp/image/'+i+'.png')
	else:
		pass
		
	imgs_tate = []
	for i in gbs_select:
		imgs_tate.append(Image.open('tmp/image/'+i+'.png'))
		
	imgs_tate_r = []
	all_width = [np.array(i).shape[0] for i in imgs_tate]
	max_len = max(all_width)
	for i in range(len(all_width)):
		sa = max_len-all_width[i]
		imgs_tate_r.append(np.array(add_margin(imgs_tate[i],sa,all_width[i])))
		
	naravi = 16
	_,amari = divmod(len(imgs_tate_r),naravi)
		
	if amari != 0:
	    imgs_tate_r = imgs_tate_r + [np.full(imgs_tate_r[0].shape,255)] * (naravi-amari)
	else:
		pass
		
	imgs_yoko = []
	for i in range(0,len(imgs_tate_r),naravi):
		imgs_yoko.append(np.concatenate(imgs_tate_r[i:i+naravi],axis=1))
		
	#縦のサイズ調整
	imgs_yoko_resized = []
	for i in imgs_yoko:
	    tn = np.where(i==255, np.nan, i) # 255をnanにする
	    not_white = []
	    for j in range(len(tn[:])):
	        if np.prod([np.isnan(tn[j,:,x]) for x in range(0,4)]) == 0: # 4次元で255だらけのところを外す
	            not_white.append(j)
	    space = 80
	    iro_ari = list(range(not_white[0]-space,not_white[0]))+not_white+list(range(not_white[-1],not_white[-1]+space))
	    imgs_yoko_resized.append(i[iro_ari, :,:])
		
	fig, ax = plt.subplots()
	ax.axis('off')
	
	im = np.concatenate(imgs_yoko_resized,axis=0)
	yoko_space = []
	for i in range(im.shape[1]-40):
	    if np.sum(im[:,i,0]-im[:,i+40,0]) == 0:
	        pass
	    else:
	        yoko_space.append(i)
	im_2 = im[:,yoko_space,:]
	
	
	pil_im = Image.fromarray(im_2.astype(np.uint8))
	pil_im.save('cMAGS.png')  
	
def page_each(gb):
	if not os.path.isfile('tmp/image/'+gb+'.png'):
	    make_map(accession=gb,store='tmp/image/')
	else:
	    pass



parser = argparse.ArgumentParser('Option to run VicMAG')


parser.add_argument('--dir',help='path to directory containing genbank files',required=True)
parser.add_argument('--plasflow',help='plasflow file',default='')
parser.add_argument('--checkv_qua',default='')
parser.add_argument('--checkv_pro',default='')
parser.add_argument('--genomad_p',help='genomad summary_plasmid',default='')
parser.add_argument('--genomad_v',default='')

args = parser.parse_args()

os.makedirs('tmp',exist_ok=True)

gbs_select = []
gbs = {}

if os.path.isdir(args.dir):
	uploaded_files = os.listdir(args.dir)
	if len(uploaded_files)>0:
		for uploaded_file in uploaded_files:
			record = SeqIO.read(args.dir+'/'+uploaded_file,'genbank')
			gbs[record.name] = record
		gbs_select = gbs_select + sort_gb_by_length(gbs)
	else:
		pass

if os.path.isfile(args.plasflow):
    pf = pd.read_table(args.plasflow,index_col=0)
    pf.to_csv('tmp/pf.csv')
    
if os.path.isfile(args.genomad_p):
    pf = pd.read_table(args.genomad_p,index_col=0)
    pf.to_csv('tmp/geno_p.csv')    
    
if os.path.isfile(args.genomad_v):
    pf = pd.read_table(args.genomad_v,index_col=0)
    pf.to_csv('tmp/geno_v.csv')  

if os.path.isfile(args.checkv_qua) and os.path.file(args.checkv_pro):
            sum_checkv = pd.read_table(args.checkv_qua,index_col=0)
            sum_checkv.to_csv('tmp/sum_checkv.csv')
            
            fasta_checkv = SeqIO.parse(args.checkv_pro,'fasta')
            SeqIO.write(fasta_checkv,'tmp/fasta_checkv.fasta','fasta')

if len(gbs_select) > 1:
#	if selected_page == gbs_select[0]:
		page_main()
#	else:
#		page_each(selected_page)
else:
	pass


