#!/opt/python/bin/python3
# coding: utf-8

import argparse
parser = argparse.ArgumentParser(description='The program takes GTF as input and finds NMD-causing local splicing events \
                                 in the form of a list of exon-exon junctions\
                                 that discriminate between an NMD transcript and coding transcripts of the same gene.\
                                 It optionally calculates PSI values of found NMD-causing local splicing events\
                                 based on split-read counts (obtained e.g by pyIPSA package)')
parser.add_argument('-g','--gtf', metavar='infile.gtf', type=str, required=True,
                    help='input GTF file')
parser.add_argument('-o','--orf', action='store_true', 
                    help='whether to annotate ORFs. Requires --ann and --genome argument to be set')
parser.add_argument('-a','--ann', metavar='annotation.gtf', type=str,
                    help='GTF file with annotated genes, transcripts and ORFs')
parser.add_argument('-G','--genome', metavar='genome.fa', type=str,
                    help='single multifasta genome file. Chromosome names must match with GTFs')
parser.add_argument('-n','--nmd', action='store_true', 
                    help='whether to predict targets of NMD in input.gtf')
parser.add_argument('-N','--nmdann', action='store_true', 
                    help='whether to predict targets of NMD in annotation.gtf')
parser.add_argument('-r','--ref', metavar='transcripts.txt or attr:value', type=str,
                    help='Either a file with reference transcript ids, each id on a new line, or a string <attribute>:<value>. In the second case transcript is chosen as reference if the value of its <attribute> contains <value>')
parser.add_argument('-p','--prefix', metavar='path/to/output/', type=str, default='./nmdj/',
                    help='Prefix for NMDj output files')
parser.add_argument('-q','--ipsa_files', metavar='file.txt', type=str,
                    help='File with paths to ipsa files containing counts of RNA-Seq split-reads aligned to junctions')
parser.add_argument('--threads', type=int,
                    help='number of threads for parallel processing')
parser.add_argument('--no_clustering', action='store_true',
                    help='skip clustering events')



#command = '''--gtf /gss/home/l.zavileisky.ext/nmdj_paper/data/simulated/unmodified/novel/SRR1079455_stringtie_novel.gtf 
#--ann /gss/home/l.zavileisky.ext/TCGA/GRCh38_dna/annotation/Homo_sapiens.GRCh38.108.chr.gtf 
#--orf --nmd --ref transcript.txt -G /gss/dplab/genomes/GRCh38/GRCh38.primary_assembly.genome.fa -p ./test_run/ 
#-q ipsa_file.txt --threads 8'''

args = parser.parse_args()

import pandas as pd
from gtfparse import read_gtf

import warnings
warnings.filterwarnings("ignore")
import subprocess
import time
import re

import nmdj_functions as nmdj
from importlib import reload
reload(nmdj)

start_time_all = time.time()

#1 load gtf
print('Loading GTF...')
if args.orf:
    if args.ann is None:
        raise Exception('--ann argument is not set, cannot annotate ORFs')
    if args.genome is None:
        raise Exception('--genome argument is not set, cannot annotate ORFs')
        
#    gtf = nmdj.load_annotation(args.gtf)
#    ann = nmdj.load_annotation(args.ann)
    gtf = read_gtf(args.gtf,result_type='pandas').rename(columns={'seqname':'chrn','feature':'type'})
    gtf = gtf.apply(pd.to_numeric, errors='ignore')
    ann = read_gtf(args.ann,result_type='pandas').rename(columns={'seqname':'chrn','feature':'type'})
    ann = ann.apply(pd.to_numeric, errors='ignore')
    
    bigann = pd.concat([gtf,ann])
else:
#    gtf = nmdj.load_annotation(args.gtf)
    gtf = read_gtf(args.gtf,result_type='pandas').rename(columns={'seqname':'chrn','feature':'type'})
    gtf = gtf.apply(pd.to_numeric, errors='ignore')
    bigann = gtf.copy()
print(f'input GTF loaded, time: {round((time.time()-start_time_all)/60,2)} min')

#!!!! remove
#bigann.chrn = 'chr'+bigann.chrn.astype('str')
#gtf.chrn = 'chr'+gtf.chrn.astype('str')


outdir = args.prefix.split('/')
outdir = '/'.join(outdir[:-1])
if outdir!='':
    subprocess.run(f"mkdir -p {outdir}", shell=True, check=True)
    print(f"mkdir -p {outdir}")


problems = []


if args.orf:
    
#assign novel transcripts to genes
    mapd = nmdj.noveltr2gene(gtf,ann)
    gtf['gene_id'] = gtf.transcript_id.map(mapd)
    p = gtf[(gtf.gene_id.isna())&(gtf.type=='transcript')][['transcript_id']]
    p.columns = ['item_id']
    p['item_type'] = 'novel_transcript'
    p['reason'] = 'unable to assign to annotated genes'
    problems.append(p.copy())
    gtf = gtf[gtf.gene_id.notna()]
    print(f'Novel transcripts assigned to genes, time: {round((time.time()-start_time_all)/60,2)} min')

#find ORFs
    orf = nmdj.orfann(gtf,ann,args.genome)
    p = gtf[(gtf.type=='transcript')&(~gtf.transcript_id.isin(orf.transcript_id.unique()))][['transcript_id']]
    p.columns = ['item_id']
    p['item_type'] = 'novel_transcript'
    p['reason'] = 'cannot find ORF'
    problems.append(p.copy())

    out = f'{args.prefix}novel_ORF.gtf'
    nmdj.write_annotation(orf,out)
    print(f'ORFs found and written to {out}, time: {round((time.time()-start_time_all)/60,2)} min')

    if not args.nmd:
        warnings.warn("--nmd is not set, although novel transcripts are given. Gonna find NDMTs anyway)")
    if args.nmdann:
        gtf = pd.concat([ann,gtf,orf])
        gtf = nmdj.process_annotation(gtf)
        gtf = nmdj.find_nmd(gtf)
    else:
        novel = pd.concat([gtf,orf])
        novel = nmdj.process_annotation(novel)
        novel = nmdj.find_nmd(novel)
        ann = nmdj.process_annotation(ann)
        gtf = pd.concat([novel,ann])
else:
    gtf = nmdj.process_annotation(gtf)
    if args.nmd:
        gtf = nmdj.find_nmd(gtf)

if len(gtf)==0:
    raise Exception('No fully annotated transcripts.\n\
    Check if any of <transcript, exon, start_codon, stop_codon> annotations is absent')



if args.nmd:
    out = f'{args.prefix}NMD_annotation.tsv'
    gtf[['transcript_id','transcript_biotype']].drop_duplicates().to_csv(out,sep='\t',index=False)
    print(f'NMDTs found and written to {out}, time: {round((time.time()-start_time_all)/60,2)} min')



#filter out genes without NMDTs or coding transcripts
u = set(gtf[gtf.transcript_biotype=='protein_coding'].gene_id)&set(gtf[gtf.transcript_biotype=='nonsense_mediated_decay'].gene_id)

gtf = gtf[(gtf.gene_id.isin(u))&
          (gtf.transcript_biotype.isin(['protein_coding','nonsense_mediated_decay']))]
print(f'GTF filtered, time: {round((time.time()-start_time_all)/60,2)} min')



#calculate frame
gtf = nmdj.calculate_frame(gtf)
print(f'reading frames calculated, time: {round((time.time()-start_time_all)/60,2)} min')

#!!!! remove
#import sys
#out = f'{args.prefix}allsites.tsv'
#gtf.to_csv(out,sep='\t',index=False)
#sys.exit()
#!!!!


#find CJs
print('Starting to search for CJs...')
jlist,noj,stat = nmdj.process_events(gtf,threads=args.threads,start_time=start_time_all)
print(f'CJs found, time: {round((time.time()-start_time_all)/60,2)} min')

#!!!! remove
#import sys
#out = f'{args.prefix}noj.tsv'
#noj.to_csv(out,sep='\t',index=False)
#sys.exit()
#!!!!

p = stat[stat.junction_tags.apply(lambda x: 'error' in x)][['nmd_id']]
p.columns = ['item_id']
p['item_type'] = 'NMDT'
p['reason'] = 'CJ searching error'
problems.append(p.copy())
for tag,reason in zip(['NMD_only', 'coding_only', 'no_junctions', 'altstart'],
                      ['no coding CJs', 'no NMD CJs','no CJs', 'alternative start-codon']):
    p = stat[stat.junction_tags.apply(lambda x: tag in x)][['nmd_id']]
    p.columns = ['item_id']
    p['item_type'] = 'NMDT'
    p['reason'] = reason
    problems.append(p.copy())



#cluster events
def func(x):
    return set.union(*x)

cl = ''
if not args.no_clustering:
    cl = 'clustered and '
    jlist['junction'] = jlist[['start','end','jtype']].apply(tuple,1)
    event_groups = jlist.groupby('gene_id').apply(nmdj.group_events).reset_index().drop('level_1',1)

    results = event_groups.merge(jlist)

    nmdid = results.groupby(['gene_id','event_num']).nmd_transcript_id.agg(set).apply(sorted).apply(lambda x: ', '.join(x))
    nmdid.name = 'nmd_id_list'

    results = results.groupby(['gene_id','event_num','chrn','strand',
                               'transcript_biotype','start','end','jtype']).agg({'transcript_id':func,'noj':max,
                                                                                 'interval_start':min,'interval_end':max}).reset_index()

    results = results.merge(nmdid.reset_index(),how='left')
    results['nmd_id_list'] = results['nmd_id_list'].fillna('')

    results['event_id'] = results.gene_id+'_'+results.event_num.astype('str')
    results.drop('event_num',axis=1,inplace=True)
else:
    results = jlist.copy()
    results['nmd_id_list'] = results.nmd_transcript_id
    results.rename(columns={'nmd_transcript_id':'event_id'},inplace=True)

    
cid = results[results.transcript_biotype=='protein_coding']
cid = cid.groupby('event_id').transcript_id.agg(func).apply(sorted).apply(lambda x: ', '.join(x))
cid.name = 'coding_id_list'
results = results.merge(cid.reset_index(),how='left')
results['coding_id_list'] = results['coding_id_list'].fillna('')

results.rename(columns={'transcript_id':'tr'},inplace=True)

cols = ['event_id','transcript_biotype','chrn','start','end','strand','jtype','tr']
outdf = results[cols].copy()
outdf[['start','end']] = outdf[['start','end']].astype('int')
outdf.tr = outdf.tr.apply(sorted).apply(lambda x: ', '.join(x))
outdf.jtype = outdf.jtype.map({'j':'junction','ir':'intron_retention'})
outdf.rename(columns={'tr':'transcripts_with_junction',
                      'jtype':'junction_type','chrn':'chromosome'},inplace=True)
out = f'{args.prefix}junctions.tsv'
outdf.to_csv(out,sep='\t',index=False)
print(f'CJs {cl}written to {out}, time: {round((time.time()-start_time_all)/60,2)} min')


#classify
results = results[results.noj==0] #get rid of events with either only coding or only NMD CJs
if args.ref is not None:
    if ':' in args.ref:
        tag,val = args.ref.split(':',1)
        print(f'Reference transcript set is specified as <{tag}>=<{val}>')
        ref = bigann[(bigann.type=='transcript')&(bigann[tag].fillna('').str.contains(val))].transcript_id.unique()
    else:
        print(f'Reference transcript set is specified as a list of transcript ids in {args.ref}')
        with open(args.ref,'r') as file:
            ref = file.read().strip().split('\n')
else:
    print('Reference transcript set is not specified. Gonna use all coding transcripts as a reference')
    ref = None

trids = nmdj.get_trids(results,bigann,ref)


res = []
for i, (key, df) in enumerate(trids.groupby('event_id')):
    v = nmdj.classify_event(df)
    res.append([key,v])
clres = pd.DataFrame(res,columns=['event_id','VIP'])

clres['event_type'] = clres.VIP.apply(nmdj.convert_etype)
results = results.merge(clres)
clres['all_CJs'] = 'Yes'
trids['all_CJs'] = 'Yes'

#classify events lacking CJs but with interval
if len(noj)>0:
    noj['nmd_id_list'] = noj.nmd_transcript_id
    noj.rename(columns={'nmd_transcript_id':'event_id','transcript_id':'tr'},inplace=True)
    cid = noj[noj.transcript_biotype=='protein_coding']
    cid = cid.groupby('event_id').tr.agg(func).apply(sorted).apply(lambda x: ', '.join(x))
    cid.name = 'coding_id_list'
    noj = noj.merge(cid.reset_index(),how='left')
    noj['coding_id_list'] = noj['coding_id_list'].fillna('')

    trids2 = nmdj.get_trids(noj,bigann,ref)
    res = []
    for i, (key, df) in enumerate(trids2.groupby('event_id')):
        v = nmdj.classify_event(df)
        res.append([key,v])
    nojres = pd.DataFrame(res,columns=['event_id','VIP'])
    nojres['event_type'] = nojres.VIP.apply(nmdj.convert_etype)
    nojres['all_CJs'] = 'No'
    trids2['all_CJs'] = 'No'
    clres = pd.concat([clres,nojres])
    trids = pd.concat([trids,trids2])

events = pd.concat([results,noj])[['gene_id','event_id',
                                   'interval_start','interval_end',
                                   'nmd_id_list','coding_id_list']].drop_duplicates().merge(clres)


out = f'{args.prefix}transcript_structures.tsv'
trids.to_csv(out,sep='\t',index=False)
print(f'Transcript structures in event regions written to {out}')

out = f'{args.prefix}classification.tsv'
events.to_csv(out,sep='\t',index=False)
print(f'Event classification written to {out}, time: {round((time.time()-start_time_all)/60,2)} min')


if args.ipsa_files is not None:
    tag_list=[]
    coef = results.groupby('event_id').apply(nmdj.assign_coefficients,tag_list).reset_index(drop=True)
    coeftag = pd.DataFrame(tag_list, columns = ['event_id','coef_tags'])

    for tag in ['pseudo_blob', 'bad_ir', 'error']:
        p = coeftag[coeftag.coef_tags.apply(lambda x: tag in x)][['event_id']]
        p.columns = ['item_id']
        p['item_type'] = 'event'
        p['reason'] = 'PSI calculation error'
        problems.append(p.copy())

    for key,tmp in coef.groupby('event_id'):
        if nmdj.get_sum(tmp)==1:
            res.append([key,'event','PSI calculation error'])

    p = pd.DataFrame(res,columns = ['item_id','item_type','reason'])
    problems.append(p.copy())
    coef = coef[~coef.event_id.isin(p.item_id)]
    
    coef[['start','end']] = coef[['start','end']].astype('int')
    sites = coef[coef.jtype=='ir']
    junctions = coef[coef.jtype=='j']

    value_vars = ['start','end']
    id_vars = [col for col in coef.columns.tolist() if col not in value_vars]
    sites = sites.melt(id_vars = id_vars, value_vars = value_vars, value_name='coord')
    sites['long_id'] = sites.chrn.astype('str')+"_"+sites.coord.astype('str')+"_"+sites.strand
    sites.drop(['coord','variable','strand'],axis=1, inplace=True)

    junctions['long_id'] = junctions.chrn.astype('str')+\
                           "_"+junctions.start.astype('str')+"_"+junctions.end.astype('str')+"_"+junctions.strand

    catalog = pd.concat([junctions,sites])
    
    #!!!! remove
    #catalog.long_id = 'chr'+catalog.long_id
    
    out = f'{args.prefix}catalog.tsv'
    catalog.to_csv(out,sep='\t',index=False)
    
    with open(args.ipsa_files,'r') as file:
        files = file.read().strip('\n').split('\n')
    sdict = dict()
    jdict = dict()
    for file in files:
        sample_name = re.split('/|\.',file)[-3]
        if 'S6' in file:
            sdict[file] = sample_name
        elif 'J6' in file:
            jdict[file] = sample_name
        else:
            raise Exception(f'not a pyIPSA file type: {file}')
    sdata = nmdj.collect_counts(sdict.keys(),ncores=args.threads,start_time=start_time_all)
    jdata = nmdj.collect_counts(jdict.keys(),ncores=args.threads,start_time=start_time_all)
    sdata.columns = sdata.columns.map(sdict)
    jdata.columns = jdata.columns.map(jdict)
    data = pd.concat([jdata,sdata])
    
    psi = nmdj.get_psi(catalog,data)
    out = f'{args.prefix}PSI.tsv'
    psi.to_csv(out,sep='\t',index=False)
    print(f'Quantification finished, results written to {out}, time: {round((time.time()-start_time_all)/60,2)} min')



out = f'{args.prefix}failed.tsv'
pd.concat(problems).to_csv(out,sep='\t',index=False)
print(f'Statistic of failures collected, results written to {out}, time: {round((time.time()-start_time_all)/60,2)} min')


print(f'Finished. Total time: {round((time.time()-start_time_all)/60,2)} min')


