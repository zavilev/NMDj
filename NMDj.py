#NMDj

import argparse
parser = argparse.ArgumentParser(description='The program takes GTF as input and finds exon-exon junctions\
                                 that discriminate between an NMD transcript and coding transcripts of the same gene.\
                                 Then it also assigns a coefficient to each junction for PSI value calculation\
                                 based on split-read counts (obtained e.g by pyIPSA package)')
parser.add_argument('-g','--gtf', metavar='infile.gtf', type=str, required=True,
                    help='GTF file containing protein-coding genes')
parser.add_argument('-n','--nmd', action='store_true', 
                    help='whether to predict targets of NMD')
parser.add_argument('-b','--bed', metavar = 'outfile.bed', 
                    help='output file with selected junctions in bed format')
parser.add_argument('-j','--junctions', metavar = 'outfile.tsv', required=True,
                    help='output file with selected junctions and assigned coefficients for PSI calculation')
parser.add_argument('-s','--stat', metavar = 'stat.tsv',
                    help='Summary statistics of NMD events')
parser.add_argument('--bed_header', action='store_true', 
                    help = 'whether to add the first line controlling track appearance in genome browser to bed file')

args = parser.parse_args()

import pandas as pd
import numpy as np
import glob
from collections import defaultdict
from collections.abc import Iterable
from itertools import chain, combinations
import sys, os
import time
STDOUT = sys.stdout
import warnings
warnings.filterwarnings("ignore")

import nmdj_functions as nmdj

print('Processing annotation...')
ann = pd.read_csv(args.gtf,header=None,sep='\t',comment='#', names=["chrn","source","type","start","end","strand","desc"])
types = {'exon', 'start_codon', 'stop_codon'}
if types-set(ann['type'].unique()):
    raise ValueError(f'The third column of input gtf should contain exon, start_codon and stop_codon types\n\
                       {", ".join(types-set(ann["type"].unique()))} is missing')
ann = ann[ann['type'].isin(types)]


trs = nmdj.process_annotation(ann)
print('Done')

if not (('transcript_biotype'  in trs.columns) or args.nmd):
    raise ValueError('Input gtf does`t have NMD annotation (transcript_biotype tag). Please, use --nmd option')

if args.nmd:
    print('Annotating transcript with NMD tag...')
    trs = nmdj.find_nmd(trs)
    print('Done')

print('Calculating frame of splice sites inside CDS...')
allsites = nmdj.calculate_frame(trs)
print('Done')

print('Selecting junctions that discriminate between coding and NMD transcripts...')
results,stat = nmdj.process_events(allsites)
print('Done')

tags = list(set(nmdj.flatten(stat.junction_tags.unique())))
for tag in tags:
    stat[tag] = stat.junction_tags.apply(lambda x: tag in x).astype('int')

aa = results.loc[(results.is_bad==0)|
                 (results.transcript_biotype=='protein_coding')].groupby('nmd_transcript_id').transcript_biotype.unique()

bb = aa.apply(lambda x: ":".join(sorted(x))).map({'nonsense_mediated_decay':'NMD_only',
                                        'nonsense_mediated_decay:protein_coding':'good',
                                        'protein_coding':'coding_only'})

stat['good_junctions'] = stat.nmd_id.map(bb.to_dict()).fillna('no_junctions')
good_ids = stat[stat['good_junctions']=='good'].nmd_id
results = results[results.nmd_transcript_id.isin(good_ids)]

print('Classify NMD events...')
event_types = results.groupby('nmd_transcript_id').apply(nmdj.classify).to_dict()
stat['event_type'] = stat.nmd_id.map(event_types).fillna('bad')
print('Done')

results.rename(columns={'transcript_id':'tr'},inplace=True)

print('Assigning coefficients for PSI calculation to selected junctions...')
tag_list=[]
coefficients = results.groupby('nmd_transcript_id').apply(nmdj.assign_coefficients,tag_list)
coeftag = pd.DataFrame(tag_list, columns = ['nmd_id','coef_tags'])
print('Done')

tags = list(set(nmdj.flatten(coeftag.coef_tags.unique())))
for tag in tags:
    if type(tag)==str: #if not, it`s an error message
        coeftag[tag] = coeftag['coef_tags'].apply(lambda x: tag in x).astype('int')

bad_tags = ['pseudo_blob', 'bad_ir', 'error']
coeftag['good_coef'] = (coeftag[[t for t in tags if t in bad_tags]].apply(sum,1)==0).astype('int')

final_stat = pd.merge(stat,coeftag,how='outer',on = 'nmd_id')

coef = coefficients.reset_index(drop=True)

#mark redundant events
print('Marking redundant events...')
coef['long_id'] = coef[['start','end','jtype']].astype('str').apply(lambda x: "_".join(x),axis=1)
nmde = coef[coef.coef!=0].groupby(['gene_id','nmd_transcript_id',
                                   'interval_start','interval_end']).long_id.agg(list).apply(sorted).apply(tuple).reset_index()

nmdred = nmde.groupby('gene_id').apply(nmdj.reduce_redundancy).drop('long_id',1)
coef = coef.merge(nmdred)

nmdred = nmdred[['nmd_transcript_id','redundant','intersect']].rename(columns = {'nmd_transcript_id':'nmd_id'})
final_stat = final_stat.merge(nmdred,how='left')
print("Done")

if args.bed:
    print('Writing bed file with selected junctions...')
    nmdj.get_bed(coef[coef.redundant!=1],args.bed, header=args.bed_header)
    print('Done')
    
if 'chr' not in coef.chrn.astype('str').iloc[0]:
    coef.chrn = 'chr'+coef.chrn.astype('str')

sites = coef[coef.jtype=='ir']
junctions = coef[coef.jtype=='j']

value_vars = ['start','end']
id_vars = [col for col in coef.columns.tolist() if col not in value_vars]
sites = sites.melt(id_vars = id_vars, value_vars = value_vars, value_name='coord')
sites['long_id'] = sites.chrn+"_"+sites.coord.astype('str')+"_"+sites.strand
sites.drop(['coord','variable','strand'],axis=1, inplace=True)
sites = sites[['long_id','jtype','transcript_biotype','nmd_transcript_id','chrn',
               'simple_id','coef','is_bad','tr','mane','gene_id','gene_name',
               'interval_start','interval_end','redundant','intersect']]
#sites.tr = sites.tr.apply(lambda x: ", ".join(sorted(x)))
junctions['long_id'] = junctions.chrn+"_"+junctions.start.astype('str')+"_"+junctions.end.astype('str')+"_"+junctions.strand
junctions = junctions[['long_id','jtype','transcript_biotype','nmd_transcript_id','chrn',
                       'simple_id','coef','is_bad','tr','mane','gene_id','gene_name',
                       'interval_start','interval_end','redundant','intersect']]
final_results = pd.concat([junctions,sites]).sort_values(['chrn','simple_id',
                                                          'transcript_biotype',
                                                          'jtype','long_id']).reset_index(drop=True).drop('chrn',axis=1)

final_results.to_csv(args.junctions, sep='\t', index=False)
if args.stat:
    final_stat.to_csv(args.stat, sep='\t', index=False)

print('Finished\n')
print('Statistics:\n')

all_nmd = len(final_stat)
good_j = sum(final_stat.good_junctions=='good')
print(f"Total number of NMD events: {all_nmd}\n")
print(f"NMD events without junction problems: {good_j}\n")
if good_j<all_nmd:
    print('Junction_problems:')
    j_problems = final_stat[final_stat.good_junctions!='good'].good_junctions.value_counts()
    if 'altstart' in final_stat.columns:
        altstart = sum(final_stat.altstart)
    for name, n in zip(j_problems.index,j_problems):
        if name=='no_junctions':
            n-=altstart
        print(f'{name}:\t{n}')
    if altstart!=0:
        print(f'Alt start:\t{altstart}')
print('\nEvent types:')
etypes = final_stat[final_stat.event_type!='bad'].event_type.value_counts()
for name, n in zip(etypes.index,etypes):
    print(f'{name}:\t{n}')
    
good_coef = int(sum(final_stat.good_coef.fillna(0)))
print(f"\nNMD events without junction and site problems: {good_coef}\n")
if good_coef<good_j:
    print('Site problems:')
    coef_problems = final_stat[(final_stat.good_coef==0)].coef_tags.apply(lambda x: tuple([i for i in x if i!='shorter_blob'])).value_counts()
    for name, n in zip(coef_problems.index,coef_problems):
        print(f'{name}:\t{n}') 


