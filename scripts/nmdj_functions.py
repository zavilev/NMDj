#This is the file containing functions for NMDj analysis

#0. Load all needed libraries

import pandas as pd
import numpy as np
from collections import defaultdict
from collections.abc import Iterable
from itertools import combinations,repeat
import time
import warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO, Seq
import concurrent.futures
import io
import subprocess
import csv

#1. read and write gtf

def desc2dict(desc):
    desc = desc.split("; ")
    dic = defaultdict(list)
    for item in desc:
        item = item.split(" ")
        dic[item[0]].append(item[1].strip('";'))
    return {key:",".join(value) for key,value in dic.items()}

def load_annotation(gtfile):
    gtf = pd.read_csv(gtfile,header=None,sep='\t',
                      comment='#',names=["chrn","source","type","start","end",'score',"strand",'frame',"desc"])

    desc = gtf.apply(lambda x: desc2dict(x.desc), axis=1, result_type="expand")
    gtf = pd.concat([gtf,desc],axis=1)
    return gtf.apply(pd.to_numeric, errors='ignore')

def format_desc(row,*attrs):
    desc = [f'{attr} "{row[attr]}"' for attr in attrs if not pd.isnull(row[attr])]
    return '; '.join(desc)+';'

def write_annotation(df,file):
    gtfcols = ['chrn', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame' ,'desc']
    for col in set(gtfcols)-set(df.columns):
        df[col] = '.'
    attrs = [i for i in df.columns if i not in gtfcols]
    df['desc'] = df.apply(format_desc,axis=1,args=attrs)
    return df[gtfcols].to_csv(file,sep='\t',index=False,header=False,quoting=csv.QUOTE_NONE)
    
#2 Assign novel transcripts to annotated genes and find ORF

def noveltr2gene(novel,annot):
    nexons = novel[novel.type=='exon']
    noveltrs = novel[novel.type=='transcript']
    anexons = annot[annot.type=='exon']
    antrs = annot[annot.type=='transcript']
    
    nexons['tclass'] = 'novel'
    anexons['tclass'] = 'annot'
    exons = pd.concat([nexons,anexons]).sort_values(['transcript_id','start','end']).reset_index(drop=True)
    
    #get junctions of all transcripts
    jstarts = exons.groupby(['tclass',
                             'gene_id',
                             'transcript_id'],sort=False).end.apply(lambda x:[np.nan]+x.tolist()[:-1])
    exons = exons.drop('end',1).rename(columns={'start':'end'})
    exons['start'] = list(flatten(jstarts))
    
    junctions = exons[exons.start.notna()][['tclass','gene_id',
                                             'transcript_id','chrn','start','end','strand']]
    junctions.start = junctions.start.astype('int')
    junctions['junction'] = junctions[['start','end']].apply(tuple,1)
    
    anngenes = junctions[junctions.tclass=='annot'].groupby('gene_id').junction.agg(set).reset_index()
    anngenes = antrs.groupby(['gene_id',
                              'chrn','strand']).agg({'start':min,'end':max}).reset_index().merge(anngenes,
                                                                                                 how='left')
    
    noveltrj = junctions[junctions.tclass=='novel'].groupby('transcript_id').junction.agg(set).reset_index()
    noveltrs = noveltrs.groupby(['transcript_id',
                                 'chrn','strand']).agg({'start':min,
                                                        'end':max}).reset_index().merge(noveltrj,how='left')
    
    noveltrs.junction = noveltrs.junction.fillna('').apply(set)
    anngenes.junction = anngenes.junction.fillna('').apply(set)
    anngenes.rename(columns={'start':'gstart','end':'gend','junction':'gjunction'},inplace=True)

    #merge results in huge table
    m = noveltrs.merge(anngenes,how='inner')
    m['l'] = np.maximum(m.start,m.gstart)
    m['r'] = np.minimum(m.end,m.gend)
    m['inter'] = m.r-m.l
    
    #filter to select parent genes for novel transcripts
    m = m[m.inter>0]
    m['gene'] = m.inter/(m.gend-m.gstart)
    m['tr'] = m.inter/(m.end-m.start)
    m['commonj'] = m.apply(lambda x: len(x.junction&x.gjunction),1)

    m = m.sort_values(['commonj','tr','gene'],ascending=False).drop_duplicates('transcript_id')
    m = m[(m.tr>0.5)]

    return m.set_index('transcript_id').gene_id.to_dict()

def orfann(novel,startann,genome):
    
#firstly get all unique annotated start codon positions
    startann = startann[startann.type=='start_codon']
    #to process start codons split by junctions
    starts = startann.groupby(['gene_id','transcript_id']).agg({'start':min,'end':max}).reset_index()
    starts = starts.drop_duplicates(['start','gene_id'])
    starts.rename(columns={'transcript_id':'start_id',
                           'start':'start_start','end':'start_end'},inplace=True)
    
#remove gene versions
    tmp = novel.copy()
    tmp.gene_id = tmp.gene_id.apply(lambda x: x.split('.')[0] if isinstance(x,str) else x)
    starts.gene_id = starts.gene_id.apply(lambda x: x.split('.')[0] if isinstance(x,str) else x)

#select starts that are in exonic regions of novel transcripts
    intls = tmp[tmp.type=='exon'].merge(starts,how='inner',on='gene_id')
    intls['l'] = (intls.start<=intls.start_end)&(intls.end>=intls.start_end)
    intls['r'] = (intls.start<=intls.start_end)&(intls.end>=intls.start_end)
    possible = intls.groupby(['transcript_id','start_id'])[['l','r']].agg(any).apply(all,1)
    intls = possible[possible].reset_index().merge(intls).drop([0,'l','r'],axis=1)

#drop 5UTR
    forward = intls[intls.strand=='+']
    reverse = intls[intls.strand=='-']

    forward = forward[forward.end>forward.start_start]
    forward.loc[forward.start<forward.start_start,'start'] = forward.loc[forward.start<forward.start_start,'start_start']
    forward.sort_values(['transcript_id','start_id','start'],ascending=True,inplace=True)

    reverse = reverse[reverse.start<reverse.start_end]
    reverse.loc[reverse.end>reverse.start_end,'end'] = reverse.loc[reverse.end>reverse.start_end,'start_end']
    reverse.sort_values(['transcript_id','start_id','start'],ascending=False,inplace=True)
    
    intls = pd.concat([forward,reverse]).reset_index(drop=True)
    
#read genome sequence
    with open(genome, "r") as genome:
        chrs = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
        
#get cds
    def get_seq(row,chrs=chrs):
        return chrs[row.chrn].seq[row.start-1:row.end]
    
    intls['sequence'] = intls.apply(get_seq,1)
    intls.loc[intls.strand=='-','sequence'] = intls.loc[intls.strand=='-','sequence'].apply(lambda x: x.reverse_complement())
    intls.sequence = intls.sequence.apply(str)
    cds = intls.groupby(['transcript_id','start_id']).sequence.agg(lambda x: "".join(x)).apply(lambda x: Seq.Seq(x).translate())
    cds = cds.apply(str).reset_index()
    cds = intls[['transcript_id','strand']].drop_duplicates().merge(cds)
    cds['ntfromstart'] = cds.sequence.apply(lambda x: len(x.split('*')[0]))*3+3
#leave only the longest cds
    cds = cds.sort_values('ntfromstart',ascending=False).drop_duplicates('transcript_id')
    intls = cds[['transcript_id','start_id','ntfromstart']].merge(intls)

#find stop codon positions
    intls[['start','end']] = intls[['start','end']].astype('int')
    intls['length'] = intls['end'] - intls['start'] + 1
    intls.length = intls.groupby('transcript_id').length.transform(np.cumsum)

    intls['delta'] = intls.length-intls.ntfromstart
    stops = intls[intls.delta>0].groupby('transcript_id').head(1).reset_index(drop=True).copy()

    stops['a'] = stops.end-stops.delta
    stops['b'] = stops.start+stops.delta
    stops.loc[stops.strand=='-','stop_start'] = stops.loc[stops.strand=='-','b']
    stops.loc[stops.strand=='+','stop_end'] = stops.loc[stops.strand=='+','a']

    stops.loc[stops.strand=='-','stop_end'] = stops.loc[stops.strand=='-','stop_start']+2
    stops.loc[stops.strand=='+','stop_start'] = stops.loc[stops.strand=='+','stop_end']-2

    stops[['stop_start','stop_end']] = stops[['stop_start','stop_end']].astype('int')
    stops.drop(columns=['a','b'],inplace=True)

#unite starts and stops in single table
    stops = stops[['chrn', 'source', 'type', 
                   'strand', 'gene_id', 'transcript_id', 
                   'stop_start', 'stop_end']].rename(columns={'stop_start':'start', 'stop_end':'end'})
    stops.type = 'stop_codon'

    starts = startann[(startann.type=='start_codon')].rename(columns={'transcript_id':
                                                                      'start_id'}).merge(cds[['start_id',
                                                                                              'transcript_id']])

    orf = pd.concat([starts[stops.columns],stops])
    orf['source'] = 'NMDj'
    return orf

#3 convert annotation to more convenient format

def process_annotation(an):  
    all_transcripts = an[an.type.isin(['transcript','exon','start_codon','stop_codon'])].groupby('type').transcript_id.agg(set)
    u1 = set.intersection(*all_transcripts)
    an = an[an.transcript_id.isin(u1)]
    
#    trs = an[(an.type=='transcript')]
#    u2 = set(trs[trs.transcript_biotype=='protein_coding'].gene_id)&set(trs[trs.transcript_biotype=='nonsense_mediated_decay'].gene_id)
#    
#    an = an[(an.transcript_id.isin(u1))&(an.gene_id.isin(u2))&
#            (an.transcript_biotype.isin(['protein_coding','nonsense_mediated_decay']))]
    
    #get tables of exons, starts and stops
    trs = an[an.type=='exon'].copy()
    if 'exon_number' not in an.columns:
        trs['exon_number'] = 1
        exf = trs[(an.strand=='+')].sort_values(['transcript_id','start','end'])
        exr = trs[(an.strand=='-')].sort_values(['transcript_id','start','end'],ascending=False)
        trs = pd.concat([exf,exr])
        trs.exon_number = trs.groupby('transcript_id').exon_number.apply(np.cumsum)

#    trs = trs[['chrn', 'start', 'end', 'strand', 
#                               'gene_id','transcript_id', 'exon_number',
#                               'gene_name', 'transcript_biotype', 'mane']]
    
    borders = an[an.type.isin(['start_codon','stop_codon'])][['start', 'end', 'type',
                                                                  'strand', 'gene_id', 'transcript_id']].copy()
    borders = borders.groupby(['type','strand',
                               'gene_id','transcript_id']).agg({"start":[min,max],"end":[min,max]}).reset_index()
    borders.columns = [i[0] if i[1]=='' else "_".join(i) for i in borders.columns]
    
    starts = borders[borders.type=='start_codon'].copy()
    starts['start_coord'] = starts.start_min
    starts.loc[starts.strand=='-','start_coord'] = starts.loc[starts.strand=='-','end_max']
    starts = starts[['gene_id','transcript_id','start_coord']]

    stops = borders[borders.type=='stop_codon'].copy()
    stops['stop_coord'] = stops.end_max
    stops.loc[stops.strand=='-','stop_coord'] = stops.loc[stops.strand=='-','start_min']
    stops = stops[['gene_id','transcript_id','stop_coord']]
        
    #add info about start and stop to each transcript
    trs = trs.merge(starts, how='inner').merge(stops, how='inner')
    
    # 5utr, cds or 3utr
    trs['exon_type'] = 'cds'
    trs.loc[(trs.strand=='+')&(trs.end<trs.start_coord),'exon_type'] = '5utr'
    trs.loc[(trs.strand=='+')&(trs.start>trs.stop_coord),'exon_type'] = '3utr'
    
    trs.loc[(trs.strand=='-')&(trs.end<trs.stop_coord),'exon_type'] = '3utr'
    trs.loc[(trs.strand=='-')&(trs.start>trs.start_coord),'exon_type'] = '5utr'
    return trs

#4. Find NMD targets by 50nt rule

def find_nmd(trs,threshold=51):
    nmd = trs.copy()
    nmdf = nmd[(nmd.strand=='+')&(nmd.end>nmd.stop_coord)]
    nmdr = nmd[(nmd.strand=='-')&(nmd.start<nmd.stop_coord)]
    #select 3utr part of the transcript
    nmdf.loc[(nmdf.start<nmdf.stop_coord),'start'] = nmdf.loc[(nmdf.start<nmdf.stop_coord),'stop_coord']+1
    nmdr.loc[(nmdr.end>nmdr.stop_coord),'end'] = nmdr.loc[(nmdr.end>nmdr.stop_coord),'stop_coord']-1
    nmd = pd.concat([nmdf,nmdr])

    nmd['length'] = nmd.end-nmd.start+1
    nmd = nmd.sort_values(['transcript_id','exon_number']).reset_index(drop=True)
    nmd['cumlen'] = nmd.groupby('transcript_id').length.transform(np.cumsum)
    lastexons = nmd.groupby('transcript_id').tail(1).index
    custom_nmd = nmd.loc[(~nmd.index.isin(lastexons))&(nmd.cumlen>threshold)].transcript_id.unique()
    trs['transcript_biotype'] = trs.transcript_id.isin(custom_nmd).map({True:'nonsense_mediated_decay',
                                                                        False: 'protein_coding'})
    return trs

#5. Assign frame to each coding site

def calculate_frame(trs, filtering=True):
    cols = ['chrn', 'start', 'end', 'strand', 'gene_id',
            'transcript_id', 'exon_number', 'start_coord',
            'stop_coord', 'exon_type', 'transcript_biotype']
    trs = trs[cols]
    cds = trs[trs['exon_type'] == 'cds'].copy()

    #change first/last exon start/end to start/stop of translation
    cds.loc[(cds.strand=='+')&(cds.start<cds.start_coord),'start'] = cds.loc[(cds.strand=='+')&(cds.start<cds.start_coord),'start_coord']
    cds.loc[(cds.strand=='+')&(cds.end>cds.stop_coord),'end'] = cds.loc[(cds.strand=='+')&(cds.end>cds.stop_coord),'stop_coord']

    cds.loc[(cds.strand=='-')&(cds.start<cds.stop_coord),'start'] = cds.loc[(cds.strand=='-')&(cds.start<cds.stop_coord),'stop_coord']
    cds.loc[(cds.strand=='-')&(cds.end>cds.start_coord),'end'] = cds.loc[(cds.strand=='-')&(cds.end>cds.start_coord),'start_coord']

    cds['length'] = cds.end-cds.start+1
    cds = cds.sort_values(['transcript_id','exon_number']).reset_index(drop=True)

    cds['cumlen3'] = cds.groupby('transcript_id').length.transform(np.cumsum)
    cds['cumlen5'] = cds['cumlen3'] - cds.length
    cds[3] = cds['cumlen3']%3
    cds[5] = cds['cumlen5']%3

    aa = cds.groupby('transcript_id').tail(1)
    cds = cds[cds.transcript_id.isin(aa[aa[3]==0].transcript_id)]

    cds = cds.melt(id_vars=[i for i in cds.columns if i not in [3,5]],
                        var_name='site', value_name='frame')
    cds['coord'] = cds['start']
    cds.loc[(cds.strand=="+")==(cds.site==3),'coord'] = cds.loc[(cds.strand=="+")==(cds.site==3),'end']


    trs[5] = trs['start']
    trs.loc[trs.strand=="-",5] = trs.loc[trs.strand=="-",'end']
    trs[3] = trs['end']
    trs.loc[trs.strand=="-",3] = trs.loc[trs.strand=="-",'start']

    trs.drop('exon_type',axis=1, inplace=True)
    trs = trs.melt(id_vars=[i for i in trs.columns if i not in [3,5]], var_name='site', value_name='coord')

    trs['site_type'] = 'cds'

    trs.loc[(trs.strand=='+')&(trs.coord<trs.start_coord),'site_type'] = '5utr'
    trs.loc[(trs.strand=='+')&(trs.coord>trs.stop_coord),'site_type'] = '3utr'

    trs.loc[(trs.strand=='-')&(trs.coord<trs.stop_coord),'site_type'] = '3utr'
    trs.loc[(trs.strand=='-')&(trs.coord>trs.start_coord),'site_type'] = '5utr'

    allsites = trs.drop(['start','end'],axis=1).merge(cds.drop(['cumlen3','cumlen5',
                                                                'exon_type','length','start','end'], axis=1), how='left')
    bad_transcripts = allsites[(allsites.frame.isna())&(allsites.site_type=='cds')].transcript_id.unique()
    allsites = allsites[~allsites.transcript_id.isin(bad_transcripts)]
    if filtering:
        nmd_genes = set(allsites[allsites.transcript_biotype=='nonsense_mediated_decay'].gene_id)
        coding_genes = set(allsites[allsites.transcript_biotype=='protein_coding'].gene_id)
        allsites = allsites[allsites.gene_id.isin(nmd_genes&coding_genes)]
    allsites = allsites.sort_values(['gene_id','transcript_id','exon_number','site'],ascending=[True]*3+[False])
    allsites['junction_number'] = allsites['exon_number']
    allsites.loc[(allsites.site==5),'junction_number'] = allsites.loc[(allsites.site==5),'junction_number']-1
    last3 = allsites.groupby('transcript_id').tail(1).index
    allsites.loc[last3,'junction_number'] = 0
    return allsites

#6. Select junctions

def select_interval_ptc(ntr,ctrs,st,tags):
    #tags.append('cds')
    if st:
        func=max
        add5 = 1
        stop = ntr[ntr.coord>ntr.stop_coord].coord.iloc[0] #the end of exon containing PTC - right border
    else:
        func=min
        add5 = -1
        stop = ntr[ntr.coord<ntr.stop_coord].coord.iloc[0] #the end of exon containing PTC - right border   
    if any(ctrs.coord==stop): #this means we`re not interested in junctions containing this exon end
        stop = ntr.stop_coord.iloc[0] +add5*3 #now right border is PTC (add5 is needed for stop codons splitted by junctions)

    #rename columns of nmd and coding dfs to merge them    
    ctrs.drop(['transcript_biotype'], inplace=True, axis=1)
    ntr.drop(['transcript_biotype'], inplace=True, axis=1)
    cols = ['chrn','coord','strand','gene_id','site']
    ctrs.columns = [i if i in cols else 'c_'+i for i in ctrs.columns]
    ntr.columns = [i if i in cols else 'n_'+i for i in ntr.columns]
    #merge nmd with coding to find the last common frame of nmd and any coding transcript
    m = ctrs.merge(ntr,how='outer')
    temp = m[(m.n_frame==m.c_frame)&(m.c_site_type=='cds')&(m.n_site_type=='cds')]
    if len(temp)==0: #we didnt find any site with common frame, need to check start
        if any(ctrs.c_start_coord==ntr.n_start_coord.iloc[0]):
            lasteq = ntr.n_start_coord.iloc[0]
            tags.append('fromstart')
            return lasteq,stop
        else:
            #print(f"alternative start with frameshift (or no common site till stop) in {ntr.n_transcript_id.iloc[0]}")
            #baddict['altstart'].append(ntr.n_transcript_id.iloc[0])
            tags.append('altstart')
            return None
    else:
        lasteq = func(temp.coord) #here it is, the last site with common frame
        if temp[temp.coord==lasteq].site.iloc[0]==5: #remove this point from the interval where we will search for junctions
            lasteq+=add5
        return lasteq,stop

    return lasteq,stop #last coordinate with equal frame 


def select_interval_utr(ntr,ctrs,st,tags):
    tags.append('3utr')
    #interval is (stop, end of the shortest 3UTR of coding transcripts with the same stop)
    ptc = ntr.stop_coord.iloc[0]
    ends = ctrs[ctrs.stop_coord==ptc].groupby('transcript_id').tail(1).coord.tolist()
    if st:
        return ptc,min(ends)
    else:
        return ptc,max(ends)

def get_junctions(ntr,ctrs,st,tags):
    
    #get borders of the interval where we gonna search for causal junctions
    if any(ctrs.stop_coord==ntr.stop_coord.iloc[0]): #it`s not PTC, need to look at 3UTR
        #print('3utr')
        left,right =  select_interval_utr(ntr,ctrs,st,tags)
        #maybe frameshift is also important to trigger NMD
        before_stop = select_interval_ptc(ntr.copy(),ctrs.copy(),st,tags)
        if (before_stop is not None) and ('fromstart' not in tags): #there exists last common junction with common frame in NMD and coding:
            left = before_stop[0]
    else:
        aa = select_interval_ptc(ntr.copy(),ctrs.copy(),st,tags)
        if aa is None:
            return None
        else:
            left,right = aa
    #print(left,right)
    #transform table to junctions
    aa = ctrs[ctrs.junction_number!=0].pivot_table(index=['transcript_id','junction_number'],
                                                   columns='site',values='coord').reset_index().sort_values(['transcript_id',                                                                                                             'junction_number'])
    ctjs = ctrs[[i for i in ctrs.columns if i not in ['coord','site',
                                                    'exon_number']]].merge(aa,how='inner').drop_duplicates(['transcript_id',
                                                                                                            'junction_number'])

    bb = ntr[ntr.junction_number!=0].pivot_table(index=['transcript_id','junction_number'],
                                                   columns='site',values='coord').reset_index().sort_values(['transcript_id',                                                                                                         'junction_number'])
    ntj = ntr[[i for i in ntr.columns if i not in ['coord','site',
                                                    'exon_number']]].merge(bb,how='inner').drop_duplicates(['transcript_id',
                                                                                                            'junction_number'])
    #select junctions that intersect with interval
    if len(ctjs)!=0:        
        if st:
            ctjs = ctjs[(ctjs[5]>=left)&(ctjs[3]<=right)]
        else: 
            ctjs = ctjs[(ctjs[5]<=left)&(ctjs[3]>=right)]
        #drop duplicated junctions
        ctjs = ctjs.groupby(['chrn','strand','gene_id','transcript_biotype',3,5]).transcript_id.agg(set).reset_index()
    if len(ntj)!=0:
        if st:
            ntj = ntj[(ntj[5]>=left)&(ntj[3]<=right)]
        else: 
            ntj = ntj[(ntj[5]<=left)&(ntj[3]>=right)]
        #ntj = ntj[ctjs.columns]
        ntj.transcript_id = ntj.transcript_id.apply(lambda x: {x})
    
    if len(ntj)*len(ctjs)!=0:
        #drop junctions present in both coding and nmd isoforms. Will do it in external function
        res = pd.concat([ctjs,ntj])#.drop_duplicates([3,5],keep=False)
    elif len(ntj)==0:
        tags.append('likelyIR')
        res = ctjs
    elif len(ctjs)==0:
        tags.append('likelyIR')
        res=ntj
    else:
        return None
    res['nmd_transcript_id'] = ntr.transcript_id.iloc[0] #it`s unique id of nmd event
    #ntj = ntj[ctjs.columns]
    if st:
        res.rename(columns={3:'start',5:'end'}, inplace=True)
    else:
        res.rename(columns={5:'start',3:'end'}, inplace=True)
    res['jtype'] = 'j'
    res['interval_start'] = min(left,right)
    res['interval_end'] = max(left,right)
    if len(res)==0:
        res['start'] = np.nan
        res['end'] = np.nan
    return res[['chrn', 'strand', 'gene_id', 'transcript_biotype', 'start',
                'end', 'transcript_id', 'nmd_transcript_id', 'jtype', 'interval_start','interval_end']]

#finds intron retention
def get_ir(ntr,ctrs,res,st,tags): #res is output of get_junctions
    cols = ['chrn','strand','gene_id','transcript_biotype',
            'start','end','transcript_id','nmd_transcript_id','jtype','interval_start','interval_end']
    exons = pd.concat([ntr,ctrs]).pivot(index=['gene_id','transcript_id', 'transcript_biotype', 'exon_number'],
                                        columns='site',values='coord').reset_index().drop_duplicates(['transcript_id',
                                                                                                      'transcript_biotype',3,5])
    if st:
        exons.rename(columns={3:'exonend',5:'exonstart'}, inplace=True)
        #res.rename(columns={3:'jstart',5:'jend','transcript_id':'transcript_id_set'}, inplace=True)
    else:
        exons.rename(columns={3:'exonstart',5:'exonend'}, inplace=True)
        #res.rename(columns={3:'jend',5:'jstart','transcript_id':'transcript_id_set'}, inplace=True)
    res.drop(['transcript_biotype','transcript_id'],inplace=True,axis=1) 
    ir = exons.merge(res,how='outer',on='gene_id')
    ir['jtype'] = 'ir'
    ir = ir[(ir.exonstart<ir.start)&(ir.exonend>ir.end)][cols] #find exons which contain whole junctions
    #ir = ir.sample(2)[cols]
    ir = ir.groupby([col for col in cols if col != 'transcript_id']).transcript_id.agg(set).reset_index()
    ir = ir.drop_duplicates(['start','end'],keep=False)
    if len(ir)!=0:
        tags.append('IR')
    return ir

def process_event(ntr,ctrs): 
    flag=True #True means that either both NMD and coding CJs present or  get_junctions returned None. 
              #False means that some CJs are absent, yet everything is OK with interval
    nmd_id = ntr.transcript_id.iloc[0]
    tags = []
    st = ntr.strand.iloc[0]=="+"
    res = get_junctions(ntr,ctrs,st,tags)
    if res is None:
         return None, flag, [nmd_id, tuple(tags)]
    junctions = res.drop_duplicates(['start','end'],keep=False)
    if len(junctions)==0:
        flag=False
        tags.append('no_junctions')
        return res, flag, [nmd_id, tuple(tags)]
    ir = get_ir(ntr,ctrs,junctions.copy(),st,tags)
    junctions = pd.concat([junctions,ir])
    tr_types = junctions.transcript_biotype.unique()
    if len(tr_types)!=2:
        flag = False
        if tr_types[0]=='nonsense_mediated_decay':
            tags.append('NMD_only')
        else:
            tags.append('coding_only')
        return res, flag, [nmd_id, tuple(tags)]
    return junctions, flag, [nmd_id, tuple(tags)]

def process_events(allsites,threads=5,start_time=None):
    results = []
    noj = []
    statistics = []
    i=0
    if start_time is None:
        start_time = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = dict()
        for gene, temp in allsites.groupby('gene_id'):
            ntrs = temp[temp.transcript_biotype=='nonsense_mediated_decay']
            ctrs = temp[temp.transcript_biotype=='protein_coding']
            for transcript,ntr in ntrs.groupby('transcript_id'):
                futures[executor.submit(process_event, ntr,ctrs)] = transcript
        print(f'tasks created, time: {round((time.time()-start_time)/60, 2)} min')
        for future in concurrent.futures.as_completed(futures):
                try:
                    res,flag,tags = future.result()
                except Exception as e:
                    print(futures[future])
                    print(e)
                    statistics.append([futures[future],('error')])
                    continue
                statistics.append(tags)
                if res is None:
                    continue
                if flag:
                    res['noj'] = 0
                    results.append(res.copy())
                else:
                    res['noj'] = 1
                    noj.append(res.copy())
                    res = res.drop_duplicates(['start','end','jtype'],keep=False)
                    results.append(res.copy())
                i+=1
                if i%1000==0:
                    print(f'Done {i} transcripts, time: {round((time.time()-start_time)/60, 2)} min')
    results = pd.concat(results)
    if len(noj)==0:
        noj = results.iloc[:0]
    else:
        noj = pd.concat(noj)
    statistics = pd.DataFrame(statistics,columns=['nmd_id','junction_tags'])
    return results,noj,statistics

#7. Classify events

#get non-redundant transcript structures in event regions
def get_trids(catalog,an,ref=None):
    cat = catalog.copy()
    if ref is not None:
        condition = catalog.tr.apply(lambda x: len(x&set(ref))>0)
        cat = cat[condition|(cat.transcript_biotype=='nonsense_mediated_decay')]
    df = cat.groupby(['event_id','nmd_id_list', 'coding_id_list']).agg({'start':min,'end':max}).reset_index()
    df.rename(columns={'start':'vmin','end':'vmax'},inplace=True)
    df.vmin-=1
    df.vmax+=1

    resdf = []
    for _,row in df.iterrows():
        for ctr in row.coding_id_list.split(', '):
            r = row.tolist()+['protein_coding',ctr]
            resdf.append(r)
        for ntr in row.nmd_id_list.split(', '):
            r = row.tolist()+['nonsense_mediated_decay',ntr]
            resdf.append(r)
    resdf = pd.DataFrame(resdf,columns=df.columns.tolist()+['transcript_biotype','transcript_id'])
    if ref is not None:
        resdf = resdf[(resdf.transcript_id.isin(ref))|(resdf.transcript_biotype=='nonsense_mediated_decay')]

    ann = resdf.merge(an[an.type=='exon'][['transcript_id','start','end','strand']])

    ann = ann[(ann.end>=ann.vmin)&(ann.start<=ann.vmax)]
    ann.loc[ann.start<=ann.vmin,'start'] = ann.loc[ann.start<=ann.vmin,'vmin']
    ann.loc[ann.end>=ann.vmax,'end'] = ann.loc[ann.end>=ann.vmax,'vmax']

    ann['exons'] = ann[['start','end']].apply(tuple,1)
    trids = ann.groupby(['event_id','strand','transcript_id','transcript_biotype']).exons.agg(list)\
               .apply(sorted).apply(tuple).reset_index()
    trids = trids.groupby(['event_id','strand','exons','transcript_biotype']).transcript_id.agg(list).reset_index()
    trids.transcript_id = trids.transcript_id.apply(sorted).apply(lambda x: ', '.join(x))
    return trids

#buld splice graph and find VIPs
class Node():
    def __init__(self,name,t=None):
        self.name = name
        self.type = t
        self.parents = []
        self.children = []
    def __str__(self):
        return f"Node {self.name} of type {self.type}\
        \nchildren: {', '.join([str(i) for i in self.children])}\
        \nparents: {', '.join([str(i) for i in self.parents])}"
    
def get_path(node1,node2,graph,strand,ind=0):
    first=True
    path = []
    curr = node1
    path.append(graph[curr].type)
    while curr!=node2:
        ind = ind if first else 0
        curr = graph[curr].children[ind]
        if curr==node1:
            return 'TASSerror'
        path.append(graph[curr].type)
        first=False
    path = ''.join(path)
    if strand=='+':
        return path
    return path[::-1]

    
def _classify(nmd,coding,strand):
    nmdi = [(nmd[i][1],nmd[i+1][0]) for i in range(len(nmd)-1)]
    codingi = [(coding[i][1],coding[i+1][0]) for i in range(len(coding)-1)]

    graph = []
    exons = sorted(set(nmd+coding))
    acc = {i[0] for i in exons}
    don = {i[1] for i in exons}
    if strand == '+':
        a = 'A'
        d = 'D'
    else:
        a = 'D'
        d = 'A'
    for n in acc:
        node = Node(n,t=a)
        graph.append(node)
    for n in don:
        node = Node(n,t=d)
        graph.append(node)
    graph.sort(key=lambda x: x.name)
    graph = {i.name:i for i in graph}

    for s,e in nmd+nmdi:
        if e not in graph[s].children:
            graph[s].children.append(e)
        if s not in graph[e].parents:
            graph[e].parents.append(s)
    for s,e in coding+codingi:
        if e not in graph[s].children:
            graph[s].children.append(e)
        if s not in graph[e].parents:
            graph[e].parents.append(s)

    left = [i for i,j in graph.items() if len(j.children)>1]
    right = [i for i,j in graph.items() if len(j.parents)>1]

    if len(left)>1:
        etype =  'complex'
    elif len(right)>1:
        etype = 'complex'
    elif (len(left)==1) and (len(right)==1):
        left = left[0]
        right = right[0]
        if left>right:
            etype =  'alt_both'
        else:
            nmdp = get_path(left,right,graph,strand,ind=0)
            cdp = get_path(left,right,graph,strand,ind=1)
            etype =  nmdp+':'+cdp
    elif (len(left)==0)==(strand=="+"):
        etype = 'alt_start'
    else:
        etype = 'alt_end'
    return etype

def classify_event(df):
    #try:
    strand = df.strand.iloc[0]
    types = set()
    for _,nmdrow in df[df.transcript_biotype=='nonsense_mediated_decay'].iterrows():
        for _,codingrow in df[df.transcript_biotype=='protein_coding'].iterrows():
            if nmdrow.exons==codingrow.exons:
                return 'same_junctions'
            nmd = list(nmdrow.exons)
            coding = list(codingrow.exons)
            etype = _classify(nmd,coding,strand)
            types|={etype}
    return '|'.join(sorted(types))
#    except Exception as e:
#        return f'error: {e}'

def _convert_etype(x):
    if x=='':
        return 'no_ref_in_region'
    if x.startswith('DA:'):
        n = int(len(x.split(':')[1])/2)-1
        if n==1:
            return 'EE'
        return f'EE{n}'
    if x.endswith(':DA'):
        n = int(len(x.split(':')[0])/2)-1
        if n==1:
            return 'PE'
        return f'PE{n}'

    if x.startswith('AD:'):
        n = int(len(x.split(':')[1])/2)-1
        if n==1:
            return 'IR'
        return f'IR{n}'
    if x.endswith(':AD'):
        n = int(len(x.split(':')[0])/2)-1
        if n==1:
            return 'ID'
        return f'ID{n}'
    human_etypes = {'DAD:DAD':'A3SS', 'ADA:ADA':'A5SS','DADA:DADA':'MXE', 'ADAD:ADAD':'A5SS+A3SS',
                    'ADA:ADADA':'A5SS+EE', 'DAD:DADAD':'EE+A3SS', 'DADAD:DAD':'PE+A3SS',
                    'ADADA:ADA':'A5SS+PE', 'ADAD:ADADAD':'A5SS+EE+A3SS', 'DAD:DADADAD':'EE2+A3SS',
                    'ADA:ADADADA':'A5SS+EE2', 'ADADADA:ADA':'A5SS+PE2', 'DAD:DADADADADAD':'EE4+A3SS',
                    'ADA:ADADADADA':'A5SS+EE3', 'DAD:DADADADAD':'EE3+A3SS', 'DAD:DADADADADADADADADADADAD':'EE10+A3SS',
                    'ADA:ADADADADADA':'A5SS+EE4'}
    if x in human_etypes.keys():
        return human_etypes[x]
    return x
def convert_etype(x):
    res = []
    for path in x.split('|'):
        res.append(_convert_etype(path))
    return '|'.join(res)

#8. Group events
        
class Event():
    def __init__(self,trid,junctions,minv,maxv):
        self.trids = [trid]
        self.junctions = junctions
        self.min = minv
        self.max = maxv
    def update(self,trid,junctions,start,end):
        self.trids.append(trid)
        self.junctions|=junctions
        self.min = min(self.min,start)
        self.max = max(self.max,end)
    def __str__(self):
        return f"Transcripts = [{', '.join(self.trids)}], min={self.min}, max={self.max}\nJunctions = {self.junctions},"
def group_events(df): #df is a set of event-junctions of the same gene 
    #df['junction'] = df[['start','end']].apply(tuple,1)
    events = dict()
    nmd = df[df.transcript_biotype=='nonsense_mediated_decay']
    for key,row in nmd.groupby('nmd_transcript_id').agg({'start':min,
                                                         'end':max,
                                                         'junction':set}).sort_values(['start','end']).iterrows():
        k=-1
        for k in events.keys():
            if (events[k].junctions & row.junction):
                events[k].update(key,row.junction,row.start,row.end)
                break
        else:
            events[k+1] = Event(key,row.junction,row.start,row.end)

    res = []
    for k,ev in events.items():
        for t in ev.trids:
            res.append([t,k])
    return pd.DataFrame(res,columns=['nmd_transcript_id','event_num'])
        
#9. Assign coefficients to selected junctions

def build_graph(df): #df is sorted by start and end
    graph = defaultdict(dict)
    #define start node
    graph['start']['tr'] = set().union(*df.tr.tolist())
    graph['start']['start'] = None
    graph['start']['end'] = min(df.start)-1
    graph['start']['parents'] = set()
    #iteratively add nodes and search for parents
    for nj,row in df.iterrows():
        trprev = {trid:(-1,'start') for trid in row.tr} #default parents
        for node,attr in graph.items():
            if (row.tr&attr['tr']) and (row.start>attr['end']): #there are common transcripts between node and added node
                for trid in row.tr&attr['tr']:
                    curr = trprev[trid][0]
                    if curr<attr['end']: #there exists rightmost junction, which can be a parent of newly added
                        trprev[trid] = (attr['end'],node) #redefine parent
#        print(f"Added node: {nj}")
#        print(trprev)
        parents = {a[1] for a in trprev.values()} #final set of parents
        for parent in parents:
            if 'children' not in graph[parent].keys(): #add node we are going to add to list of children
                graph[parent]['children'] = set()
            graph[parent]['children'].add(nj)
        #now add this node to dict of nodes
        graph[nj]['tr'] = row.tr
        graph[nj]['start'] = row.start
        graph[nj]['end'] = row.end
        graph[nj]['parents'] = parents
        graph[nj]['jtype'] = row.jtype
    #define end node
    graph['end']['tr'] = set().union(*df.tr.tolist())
    graph['end']['start'] = max(df.end)+1
    graph['end']['end'] = None
    graph['end']['children'] = set()
    graph['end']['parents'] = set()
    #link the end node to nodes that have no parents
    for node,attr in graph.items():
        if 'children' not in attr.keys():
            graph[node]['children'] = {'end'}
            graph['end']['parents'].add(node)
    return graph

    
def get_bifurcation(graph): #find blobs and sort them according to length of blob
    pairs = []
    starts = {key:attr for key,attr in graph.items() if len(attr['children'])>1}
    ends = {key:attr for key,attr in graph.items() if len(attr['parents'])>1}
    for start,sattr in starts.items():
        for end,eattr in ends.items():
            if (sattr['end']<eattr['start']):# and (len(sattr['tr']&eattr['tr'])>1):
                startset = {tuple(sorted(graph[node]['tr'])) for node in sattr['children']}
                endset = {tuple(sorted(graph[node]['tr'])) for node in eattr['parents']}
                if len(startset&endset)>1:
                    pairs.append([start,end,eattr['start']-sattr['end'],tuple(sorted(sattr['tr']&eattr['tr']))])
    pairs = sorted(pairs,key=lambda x: x[2])
    lastbif = ['start','end',graph['end']['start']-graph['start']['end'],graph['start']['tr']]
    pairs.append(lastbif)
    pairs = pd.DataFrame(pairs,columns=['start','end','length','tr']).drop_duplicates(['start','end'],
                                                                                      keep='first').reset_index(drop=True)
    delete_ind = []
    for row1,row2 in combinations(pairs.iterrows(),2):
        i1,r1 = row1
        i2,r2 = row2
        start1 = r1.start if r1.start!='start' else -1
        start2 = r2.start if r2.start!='start' else -1
        end1 = r1.end if r1.end!='end' else np.inf
        end2 = r2.end if r2.end!='end' else np.inf
        if not ((end1<=start2) or (end2<=start1)): #they intersect
            if r1.tr == r2.tr: #they belong to the same set of transcripts
                #need to delete the longer one
                if r1.length>r2.length:
                    delete_ind.append(i1)
                else:
                    delete_ind.append(i2)
    #print(delete_ind)
    pairs = pairs[~pairs.index.isin(delete_ind)] 
    pairs.tr = pairs.tr.apply(set)
    return pairs

def flatten(xs):
    for x in xs:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            yield from flatten(x)
        else:
            yield x

def get_sum(df):
    trs = list(set.union(*df.tr.tolist()))
    for tr in trs:
        s = sum(df.tr.apply(lambda x: tr in x).astype('int')*df.coef)
        if s!=1:
            return 1
    return 0

def _assign_coefficients(df): #iteratively reduce blobs, while assigning coefficients to junctions
    tags = []
    graph = build_graph(df)
    coef = defaultdict(lambda: 1)
    #get bifurcations and remove all empty paths
    empty_path=False
    while True:
        stop = True
        bifurcation = get_bifurcation(graph)
        for _,row in bifurcation.iterrows():
            if row.end in graph[row.start]['children']: #empty path
                empty_path=True
                stop=False
                #!!!!! we should remove not empty path but all the other paths in this blob
                graph[row.start]['children'] = {row.end}
                graph[row.end]['parents'] = {row.start}
                need2remove = []
                for node,attr in graph.items():
                    if (node=="start") or (node=="end"):
                        continue
                    if (attr['start']>graph[row.start]['end']) and (attr['end']<graph[row.end]['start']):
                        if attr['tr']&row.tr:
                            if attr['tr']-row.tr:
                                #it`s an error - pseudoblob
                                tags.append('pseudo_blob')
                                return None, tuple(sorted(set(tags)))
                            else:
                                need2remove.append(node)
                for node in need2remove:
                    graph.pop(node, None)
                    coef[node] = 0
        if stop:
            break
    if empty_path:
        tags.append('empty_path')

    for _,bif in bifurcation.iterrows():
        if (bif.start not in graph.keys()) or (bif.end not in graph.keys()): #this means it was reduced by another blob
            tags.append('no_bif_node')
            continue

        parallel = []
        start_paths = [child for child in graph[bif.start]['children'] if len(graph[child]['tr']&bif.tr)!=0]
        if len(start_paths)==0:
            raise Exception('no children')
        if len(start_paths)==1: #shorter blob reduced number of children
            tags.append('shorter_blob')

        for current in start_paths:#as it`s the smallest blob the number of paths = the number of children
            if len(graph[current]['tr']&bif.tr)==0:
                #it`s wrong path not leading to end
                tags.append('wrong_path')
                continue
            path = []
            curr = current
            while curr!=bif.end: #find sequential path from start to end (list of junctions)
                if len(graph[curr]['parents'])>1:
                    #print(f'wrong number of parents: {df.nmd_transcript_id.iloc[0]} {curr}')
                    tags.append('pseudo_blob')
                    return None, tuple(sorted(set(tags)))    
                path.append(curr)

                if len(graph[curr]['children'])!=1:
                    #raise Exception(f'wrong number of children: {curr}')
                    #print(f'wrong number of children: {df.nmd_transcript_id.iloc[0]} {curr}')
                    tags.append('pseudo_blob')
                    return None, tuple(sorted(set(tags)))
                last=curr
                curr = next(iter(graph[curr]['children']))

            new = tuple(flatten(path))
            parallel.append(new)
            if len(path)>1: #reduce this path to one composite junction, add coefficients
                a = 1/len(path)
                new_trs = set()
                for node in path:
                    #delete node from graph
                    new_trs |= graph[node]['tr']
                    graph.pop(node,None)
                    #add coeffcient
                    for n in flatten([node]):
                        coef[n]*=a
                #redefine children and parents of start and end nodes
                newattr = {'parents':{bif.start},'children':bif.end,'tr':new_trs}
                graph[new] = newattr
                start_children = graph[bif.start]['children']
                graph[bif.start]['children'] = (start_children-{current}) | {new}
                end_parents = graph[bif.end]['parents']
                graph[bif.end]['parents'] = (end_parents-{last}) | {new}
        #now inside blob only parallel junctions exist
        #delete parallel junctions and create new composite node
        newpar_trs = set()
        for par in parallel:
            if len(par)==1:
                newpar_trs |= graph[par[0]]['tr']
                graph.pop(par[0],None)
                graph[bif.start]['children'].remove(par[0])
                graph[bif.end]['parents'].remove(par[0])
            else:
                newpar_trs |= graph[par]['tr']
                graph.pop(par,None)
                graph[bif.start]['children'].remove(par)
                graph[bif.end]['parents'].remove(par)
        new_par = tuple(flatten(parallel))
        newpar_attr = {'children':{bif.end},"parents":{bif.start},'tr':newpar_trs}
        graph[new_par] = newpar_attr
        #redefine children and parents of start and end nodes
        graph[bif.start]['children'].add(new_par)
        graph[bif.end]['parents'].add(new_par)
        
    df['coef'] = df.index.map(coef)
    return df,tuple(sorted(set(tags)))

def check_ir(df): #sorted and index reset
    bad_ind = set()
    if sum(df.jtype=='ir')>1:
        for i1,i2 in combinations(df[df.jtype=='ir'].index,2):
            if df.loc[i1,'end']>df.loc[i2,'start']: #this IRs intersect
                if df.loc[i1,'tr']&df.loc[i2,'tr']:
                    bad_ind|={i1,i2}
    bad_ind = sorted(bad_ind)
    return bad_ind

def assign_coefficients(df,tag_list):
    try:
        df = df.sort_values(['start','end']).reset_index(drop=True)
        trs = set.union(*df.tr)
        problem = check_ir(df)
        ir_flag = False
        while len(problem)>0:
            ir_flag = True
            df = df.drop(problem,axis=0)
            problem = check_ir(df)
        trs_new = set.union(*df.tr)
        if trs_new != trs: #deletion of problematic IR leads to inability to quantify some transcripts
            tag_list.append([df.event_id.iloc[0],('fatal_ir',)])
            return None
        nmd = df[df.transcript_biotype=='nonsense_mediated_decay'].copy().reset_index(drop=True)
        coding = df[df.transcript_biotype=='protein_coding'].copy().reset_index(drop=True)

        coding,tags = _assign_coefficients(coding)
        nmd,ntags = _assign_coefficients(nmd)
        tags = ['coding:'+i for i in tags]
        ntags = ['nmd:'+i for i in ntags]
        tags+=ntags
        if ir_flag:
            tags.append('bad_ir')

        if (coding is None) or (nmd is None):
            tag_list.append([df.event_id.iloc[0],tags])
            return None

        res = pd.concat([nmd,coding])
        res.loc[res.jtype=='ir','coef'] = res.loc[res.jtype=='ir','coef']/2
        tag_list.append([df.event_id.iloc[0],tags])
        return res
    except Exception as e:
        print('error in',df.event_id.iloc[0])
        tag_list.append([df.event_id.iloc[0],('error',e)])
        return None
    
#10 parse count files and calculate PSI

def get_counts(file, se):
    #filename = file.split("/")[-1].split(".")[0]
    df = pd.read_csv(file, sep = "\t", header = None, usecols=[0,1], names=["SE",file])
    out = pd.merge(se,df, how = "left", on = "SE").fillna(0)
    out.index = out.SE
    return out.drop('SE',axis=1)

def read_shell(command, **kwargs):
    proc = subprocess.Popen(command, 
                            shell=True,
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
    output, error = proc.communicate()

    if proc.returncode == 0:
        with io.StringIO(output.decode()) as buffer:
            return pd.read_table(buffer, **kwargs)
    else:
        raise Exception(error.decode())

def collect_counts(files, ncores=None, se=None, step=10, start_time=None):
    if se is None:
        command = f'zcat {" ".join(files)} | cut -f1 | sort -u'
        se = read_shell(command, header=None, names=['SE'])
    counts = []
    if start_time is None:
        start_time = time.time()
    print("#-----starting parallel extraction of counts")
    with concurrent.futures.ProcessPoolExecutor(max_workers=ncores) as executor:
        it=0
        for out in executor.map(get_counts, files, repeat(se)):
            counts.append(out)
            if it%step==0:
                print(f"{it} files processed, time since start: {round((time.time()-start_time)/60,2)} min")
            it+=1
    counts = pd.concat(counts,axis=1)
    return counts

def get_psi(catalog,data,event_col = 'event_id', 
            nmd_col = 'transcript_biotype', 
            nmd_tag = 'nonsense_mediated_decay',sample_order=None):
    samplecols = data.columns.tolist()
    data = data.reset_index().rename(columns={'SE':'long_id'})
    data = catalog.merge(data,how='inner')
    meta = data[[i for i in data.columns if i not in samplecols]]
    arr = data[samplecols].fillna(0).mul(data.coef, axis=0)
    data = pd.concat([meta,arr],axis=1)    
    
    data[nmd_col] = (data[nmd_col]==nmd_tag).map({True: 'nmd', False: 'coding'})
    data = data.groupby([event_col,nmd_col])[arr.columns].agg(sum).reset_index()
    data = data.melt(id_vars = [event_col,nmd_col], 
                     var_name = 'sample_num', value_name='count').pivot(index=[event_col,'sample_num'], 
                                                                       columns=nmd_col).reset_index().fillna(0)
                                                                   
    data.columns = [i[0] if i[1]=='' else i[1] for i in data.columns]
    data['denom'] = data.coding+data.nmd
    data['psi'] = data.nmd/data.denom
    data.drop(['coding'], axis=1, inplace=True)
    
    if sample_order is not None:
        sample_order = pd.read_csv(sample_order)
        data = data.merge(sample_order).drop('sample_num', axis=1)
    else:
        data.rename(columns={'sample_num':'sample_id'},inplace=True)
    return data
