#This is the file containing functions for NMDj analysis

#0. Load all needed libraries

import pandas as pd
import numpy as np
from collections import defaultdict
from collections.abc import Iterable
from itertools import combinations
import sys, os
import time
STDOUT = sys.stdout
import warnings
warnings.filterwarnings("ignore")

#1. Process gtf

def desc2dict(desc):
    desc = desc.split("; ")
    dic = defaultdict(list)
    for item in desc:
        item = item.split(" ")
        dic[item[0]].append(item[1].strip('";'))
    return {key:",".join(value) for key,value in dic.items()}

#here we assume that "an" is pandas.DataFrame with columns "chrn","source","type","start","end","strand","desc"
#"chrn" is chromosome name, "desc" is semicolon-separated list of tag-value pairs
def process_annotation(an):
#    an = an[(an.type.isin(['exon',
#                            'start_codon','stop_codon']))&(an.desc.str.contains('protein_coding|nonsense_mediated_decay'))]
    
    aa = an.apply(lambda x: desc2dict(x.desc), axis=1, result_type="expand")
    an = pd.concat([an,aa],axis=1)   
    an = an[~(an.tag.str.contains('mRNA_start_NF')|an.tag.str.contains('mRNA_end_NF'))]
    all_transcripts = an.groupby('type').transcript_id.agg(set)#.tolist()
    u = set.intersection(*all_transcripts)
    an = an[(an.transcript_id.isin(u))&(an.transcript_biotype.isin(['protein_coding','nonsense_mediated_decay']))]
    an['mane'] = an.tag.fillna('').str.contains('MANE_Select').astype('int')
    an.exon_number = an.exon_number.astype('int')

    #get tables of exons, starts and stops
    trs = an[an.type=='exon'][['chrn', 'start', 'end', 'strand', 
                               'gene_id','transcript_id', 'exon_number',
                               'gene_name', 'transcript_biotype', 'mane']].copy()
    
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

#2. Find NMD targets by 50nt rule

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

#3. Assign frame to each coding site

def calculate_frame(trs, filtering=True):
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

#4. Select junctions

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
    ntr.drop(['transcript_biotype','mane'], inplace=True, axis=1)
    cols = ['chrn','coord','strand','gene_id','gene_name','mane','site']
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

def get_junctions(ntr,ctrs,st,tags, MANE):
    
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
    #return ntj,ctjs
    if len(ctjs)!=0:        
        if st:
            ctjs = ctjs[(ctjs[5]>=left)&(ctjs[3]<=right)]
        else: 
            ctjs = ctjs[(ctjs[5]<=left)&(ctjs[3]>=right)]
        #drop duplicated junctions
        ctjs = ctjs.groupby(['chrn','strand','gene_id',
                             'gene_name','transcript_biotype',3,5]).agg({'mane':max,'transcript_id':set}).reset_index()
    if len(ntj)!=0:
        if st:
            ntj = ntj[(ntj[5]>=left)&(ntj[3]<=right)]
        else: 
            ntj = ntj[(ntj[5]<=left)&(ntj[3]>=right)]
        #ntj = ntj[ctjs.columns]
        ntj.transcript_id = ntj.transcript_id.apply(lambda x: {x})
    
    if len(ntj)*len(ctjs)!=0:
        #drop junctions present in both coding and nmd isoforms
        if MANE:
            nmd_j = ntj[[3,5]].apply(tuple,axis=1)
            mane_j = ctjs[ctjs.mane==1][[3,5]].apply(tuple,axis=1)
            ctjs = ctjs[~ctjs[[3,5]].apply(tuple,axis=1).isin(nmd_j)]
            ntj = ntj[~ntj[[3,5]].apply(tuple,axis=1).isin(mane_j)]
            res = pd.concat([ntj,ctjs])
        else:    
            res = pd.concat([ntj,ctjs]).drop_duplicates([3,5],keep=False)
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
    if 'gene_name' not in res.columns:
        res['gene_name'] = pd.NA
    res['interval_start'] = min(left,right)
    res['interval_end'] = max(left,right)
    if len(res)==0:
        res['start'] = np.nan
        res['end'] = np.nan
    return res[['chrn', 'strand', 'gene_id', 'gene_name', 'transcript_biotype', 'start',
                'end', 'mane', 'transcript_id', 'nmd_transcript_id', 'jtype', 'interval_start','interval_end']]

#finds intron retention
def get_ir(ntr,ctrs,res,st,tags): #res is output of get_junctions
    cols = ['chrn','strand','gene_id','gene_name','transcript_biotype',
            'start','end','mane','transcript_id','nmd_transcript_id','jtype','interval_start','interval_end']
    exons = pd.concat([ntr,ctrs]).pivot(index=['gene_id','transcript_id', 'transcript_biotype', 'exon_number','mane'],
                                        columns='site',values='coord').reset_index().drop_duplicates(['transcript_id',
                                                                                                      'transcript_biotype',3,5])
    if st:
        exons.rename(columns={3:'exonend',5:'exonstart'}, inplace=True)
        #res.rename(columns={3:'jstart',5:'jend','transcript_id':'transcript_id_set'}, inplace=True)
    else:
        exons.rename(columns={3:'exonstart',5:'exonend'}, inplace=True)
        #res.rename(columns={3:'jend',5:'jstart','transcript_id':'transcript_id_set'}, inplace=True)
    res.drop(['transcript_biotype','transcript_id','mane'],inplace=True,axis=1) 
    ir = exons.merge(res,how='outer',on='gene_id')
    ir['jtype'] = 'ir'
    ir = ir[(ir.exonstart<ir.start)&(ir.exonend>ir.end)][cols] #find exons which contain whole junctions
    #ir = ir.sample(2)[cols]
    ir = ir.groupby([col for col in cols if col not in  ['transcript_id','mane']]).agg({'mane':max,
                                                                                        'transcript_id':set}).reset_index()
    ir = ir.drop_duplicates(['start','end'],keep=False)
    if len(ir)!=0:
        tags.append('IR')
    return ir

def process_event(ntr,ctrs, MANE):
    nmd_id = ntr.transcript_id.iloc[0]
    tags = []
    st = ntr.strand.iloc[0]=="+"
    junctions = get_junctions(ntr,ctrs,st,tags, MANE)
    if junctions is None:
         return None, [nmd_id, tuple(tags)]
    if len(junctions)==0:
        tags.append('no_junctions')
        return None, [nmd_id, tuple(tags)]
    ir = get_ir(ntr,ctrs,junctions.copy(),st,tags)
    res = pd.concat([junctions,ir])
    tr_types = res.transcript_biotype.unique()
    if len(tr_types)!=2:
        if tr_types[0]=='nonsense_mediated_decay':
            tags.append('NMD_only')
        else:
            tags.append('coding_only')
    return res, [nmd_id, tuple(tags)]

def process_events(allsites, MANE=False):
    results = []
    statistics = []
    i=0
    start_time = time.time()
    for gene, temp in allsites.groupby('gene_id'):
        ntrs = temp[temp.transcript_biotype=='nonsense_mediated_decay']
        ctrs = temp[temp.transcript_biotype=='protein_coding']
        nmdj = set(ntrs[ntrs.junction_number!=0].groupby(['transcript_id',
                                                      'junction_number']).coord.agg(list).apply(sorted).apply(tuple))
        codingj = set(ctrs[ctrs.junction_number!=0].groupby(['transcript_id',
                                                         'junction_number']).coord.agg(list).apply(sorted).apply(tuple))
        badj = nmdj&codingj
        for transcript,ntr in ntrs.groupby('transcript_id'):
            try:
                res,tags = process_event(ntr,ctrs, MANE)
            except Exception as e:
                print(transcript)
                print(e)
                statistics.append([transcript,('error')])
                continue
            statistics.append(tags)
            if res is None:
                continue
            res['simple_id'] = i
            res['is_bad'] = res[['start','end']].apply(tuple,axis=1).isin(badj).astype('int')
            results.append(res.copy())
            i+=1
            if i%1000==0:
                print(f'Done {i} transcripts, time: {round((time.time()-start_time)/60, 2)} min')
    return pd.concat(results), pd.DataFrame(statistics,columns=['nmd_id','junction_tags'])

#5. Classify events

def classify(df):
    if len(df.transcript_biotype.unique())<2:
        return 'bad'
    if len(df)>3:
        return 'complex'
    if len(df[(df.transcript_biotype=='protein_coding')&(df.mane==0)])!=0: #another junctions
        return 'complex'
    if len(df)==3: #think about cassette exon
        if any(df.jtype=='ir'):
            return 'complex'
        df = df.sort_values(['start','end']).reset_index(drop=True)
        if df.loc[0,'start']!=df.loc[1,'start']:
            return 'complex'
        if df.loc[1,'end']!=df.loc[2,'end']:
            return 'complex'
        if df.loc[1,'transcript_biotype'] == 'protein_coding': #poison exon
            return 'PE'
        else:
            return 'EE'
    if len(df)==2: #alt sites or ir
        if len(df[['start','end']].drop_duplicates())==1:
            return 'IR'
        if df.start.iloc[0] == df.start.iloc[1]:
            if df.strand.iloc[0]=="+":
                return 'A3SS'
            else:
                return 'A5SS'
        elif df.end.iloc[0] == df.end.iloc[1]:
            if df.strand.iloc[0]=="+":
                return 'A5SS'
            else:
                return 'A3SS'
        else:
            return 'complex'

#6. Get bed file of selected junctions

def get_bed(df, file, header=True):
    results = df.copy()
    results = results.sort_values(['simple_id','start','end'])
    results['score'] = 0
    results.start = results.start-1
    results['thickStart'] = results.start
    results['thickEnd'] = results.end
    results['itemRgb'] = results.transcript_biotype.map({'nonsense_mediated_decay':'0,0,0','protein_coding':'0,128,0'})
    results['blockCount'] = 2
    results['blockSizes'] = '1,1'
    results['blockStarts'] = '0,'+(results.end-results.start-1).astype('str')
    results['name'] = results.simple_id.astype('str')+":"+results.jtype
    results =  results[['chrn','start','end','name','score','strand',
                        'thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']]
    mode='w'
    if header:
        mode = 'a'
        with open(file,'w') as out:
            out.write('track name="NMD junctions" description="Junctions determining PSI" visibility=2 itemRgb=1 useScore=0\n')
    results.to_csv(file, header=False, sep='\t', index=False, mode=mode)
        
#7. Assign coefficients to selected junctions

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

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = STDOUT
    
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
                #!!!!!
                #remove empty path and then recalculate bifurcations
                #graph[row.start]['children'].remove(row.end)
                #graph[row.end]['parents'].remove(row.start)
        if stop:
            break
    if empty_path:
        tags.append('empty_path')
#    if len(bifurcation)==0:
#        df['coef'] = 1/len(df)
#        tags.append('no_bifurcation')
#        return df,tuple(sorted(tags))
    for _,bif in bifurcation.iterrows():
        if (bif.start not in graph.keys()) or (bif.end not in graph.keys()): #this means it was reduced by another blob
            tags.append('no_bif_node')
            continue
        print(f"Blob: {bif.start} - {bif.end}")
        print(f"Children: {graph[bif.start]['children']}")
        parallel = []
        start_paths = [child for child in graph[bif.start]['children'] if len(graph[child]['tr']&bif.tr)!=0]
        if len(start_paths)==0:
            raise Exception('no children')
        if len(start_paths)==1: #shorter blob reduced number of children
            print('shorter blob reduced number of children')
            print(f"Graph state: {graph.keys()}")
            tags.append('shorter_blob')
            #continue
        for current in start_paths:#as it`s the smallest blob the number of paths = the number of children
            if len(graph[current]['tr']&bif.tr)==0:
                #it`s wrong path not leading to end
                print('it`s wrong path not leading to end')
                tags.append('wrong_path')
                continue
            path = []
            curr = current
            while curr!=bif.end: #find sequential path from start to end (list of junctions)
                if len(graph[curr]['parents'])>1:
                    print(f'wrong number of parents: {df.nmd_transcript_id.iloc[0]} {curr}')
                    tags.append('pseudo_blob')
                    return None, tuple(sorted(set(tags)))    
                path.append(curr)
                #print(curr)
#                if (bif.end!='end') and (curr=='end'): #this means we walked along wrong path, which not leads to end
#                    wrongpath=True
#                    raise Exception('wrong path')
#                    break
                if len(graph[curr]['children'])!=1:
                    #raise Exception(f'wrong number of children: {curr}')
                    print(f'wrong number of children: {df.nmd_transcript_id.iloc[0]} {curr}')
                    tags.append('pseudo_blob')
                    return None, tuple(sorted(set(tags)))
                last=curr
                #curr = graph[curr]['children'].pop()
                curr = next(iter(graph[curr]['children']))
#            if wrongpath: #skip wrong path
#                continue
#            if len(path)==0:
#                #empty path
#                tags.append('empty_path')
#                #break this connection
#                graph[bif.start]['children'].remove(bif.end)
#                graph[bif.end]['parents'].remove(bif.start)
#                continue
#                print('empty path')
            new = tuple(flatten(path))
            parallel.append(new)
            print(f"Sequential: {new}")
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
        #print(graph.keys())
        #print(f"Parallel paths: {parallel}")
        #delete parallel junctions and create new composite node
        newpar_trs = set()
        for par in parallel:
            print(f"Parallel paths: {par}")
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
        #print(f"Name of node: {new_par}")
        newpar_attr = {'children':{bif.end},"parents":{bif.start},'tr':newpar_trs}
        graph[new_par] = newpar_attr
        #redefine children and parents of start and end nodes
        graph[bif.start]['children'].add(new_par)
        graph[bif.end]['parents'].add(new_par)
        print(f"Graph state: {graph.keys()}")
        for node,attr in graph.items():
            print(attr['parents'],node,attr['children'])
    df['coef'] = df.index.map(coef)
    return df,tuple(sorted(set(tags)))

def check_ir(df): #sorted and index reset
    if sum(df.jtype=='ir')>1:
        for i1,i2 in combinations(df[df.jtype=='ir'].index,2):
            if df.loc[i1,'end']>df.loc[i2,'start']: #this IRs intersect
                if df.loc[i1,'tr']&df.loc[i2,'tr']:
                    return True
    return False

def assign_coefficients(df,tag_list):
    #print(df.nmd_transcript_id.iloc[0], df.simple_id.iloc[0])
    blockPrint()
    try:
        df = df.sort_values(['start','end']).reset_index(drop=True)
        coding = df[df.transcript_biotype=='protein_coding'].copy().reset_index(drop=True)
        if check_ir(df):
            tag_list.append([df.nmd_transcript_id.iloc[0],('bad_ir',)])
            return None
        nmd = df[df.transcript_biotype=='nonsense_mediated_decay'].copy().reset_index(drop=True)
        coding,tags = _assign_coefficients(coding)
        if coding is None:
            tag_list.append([df.nmd_transcript_id.iloc[0],tags])
            return None
        nmd['coef'] = 1/len(nmd)
        res = pd.concat([nmd,coding])
        res.loc[res.jtype=='ir','coef'] = res.loc[res.jtype=='ir','coef']/2
        tag_list.append([df.nmd_transcript_id.iloc[0],tags])
        enablePrint()
        return res
    except Exception as e:
        print(df.nmd_transcript_id.iloc[0])
        tag_list.append([df.nmd_transcript_id.iloc[0],('error',e)])
        return None
    
def reduce_redundancy(df):
    intersect = set()
    red = set()
    for r1,r2 in combinations(df.itertuples(),2):
        if r1.Index == r2.Index:
            continue
        if not ((r1.interval_start>r2.interval_end) or (r2.interval_start>r1.interval_end)): #intersection
            intersect.add(r1.Index)
            intersect.add(r2.Index)
            if r1.long_id == r2.long_id: #total redundancy
                red.add(r2.Index)
    df['redundant'] = 0
    df['intersect'] = 0
    df.loc[intersect,'intersect'] = 1
    df.loc[red,'redundant'] = 1
    return df

def _classify_mane(df): #requires df with junctions of two isoforms: nmd and mane
    if len(df.transcript_biotype.unique())<2:
        return 'bad'
    if len(df[['start','end']].drop_duplicates())==1 and sum(df.jtype=='ir')==1: #ir
        return 'IR'
    if sum(df.jtype=='ir')>0:
        return 'complex'
#ASS
    if len(df)==2:
        if df.start.iloc[0] == df.start.iloc[1]:
            if df.strand.iloc[0]=="+":
                return 'A5SS'
            else:
                return 'A3SS'
        elif df.end.iloc[0] == df.end.iloc[1]:
            if df.strand.iloc[0]=="+":
                return 'A3SS'
            else:
                return 'A5SS'
#essential exons    
    if sum(df.transcript_biotype=='nonsense_mediated_decay')==1: 
        nmd = df[df.transcript_biotype=='nonsense_mediated_decay'].iloc[0]
        coding = df[df.transcript_biotype=='protein_coding'].sort_values(['start','end']).reset_index(drop=True)
        if (nmd.start==coding.iloc[0].start) and (nmd.end==coding.iloc[-1].end):
            if len(coding)==2: #single essential exon
                return 'EE'
            else: #multiple essential exons
                return f'C{len(coding)-1}EE'
        return 'complex'
#poison exons
    if sum(df.transcript_biotype=='protein_coding')==1: 
        coding = df[df.transcript_biotype=='protein_coding'].iloc[0]
        nmd = df[df.transcript_biotype=='nonsense_mediated_decay'].sort_values(['start','end']).reset_index(drop=True)
        if (coding.start==nmd.iloc[0].start) and (coding.end==nmd.iloc[-1].end):
            if len(nmd)==2: #single poison exon
                return 'PE'
            else: #multiple essential exons
                return f'C{len(nmd)-1}PE'
        return 'complex'
#MXE            
    if set(df.transcript_biotype.value_counts())=={2}:
        aa = df.groupby('transcript_biotype').agg({'start':min, 'end':max})
        if all(aa.loc['protein_coding'] == aa.loc['nonsense_mediated_decay']):
            return 'MXE'
    return 'complex'

def classify_mane(df):
    df['event_type'] = _classify_mane(df)
    return df