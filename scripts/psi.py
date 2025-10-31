from graphviz import Digraph
import nmdj_functions as nmdj
import pandas as pd
import numpy as np
import concurrent.futures
import time

from threading import Thread
import functools
from scipy.optimize import nnls



class Node():
    def __init__(self,i,trs=set(),start=None,end=None):
        self.name = i
        self.parents = set()
        self.children = set()
        self.trs = trs
        self.start = start
        self.end = end
    def add_trs(self,trs):
        self.trs|=trs
    def remove_trs(self,trs):
        self.trs-=trs
    def add_child(self,self2):
        self.children.add(self2)
        self2.parents.add(self)
    def remove_child(self,self2):
        self.children.discard(self2)
        self2.parents.discard(self)
        
    def add_parent(self,self2):
        self.parents.add(self2)
        self2.children.add(self)
    def remove_parent(self,self2):
        self.parents.discard(self2)
        self2.children.discard(self)
    def __str__(self):
        p = f"p: [{', '.join(sorted([str(i.name) for i in self.parents]))}]"
        c = f"c: [{', '.join(sorted([str(i.name) for i in self.children]))}]"
        tr = f"tr: [{', '.join(sorted([str(i) for i in self.trs]))}]"
        return f"Node {self.name}, {p}, {c}\n{tr}\nstart={self.start}, end={self.end}"
    def from_row(i,row):
        node = Node(i,trs=row.tr,start=row.start,end=row.end)
        return node
class Graph():
    def __init__(self,df=None):
        if df is None:
            alltrs = set()
            se = None
            es = None
        else:
            alltrs = set.union(*df.tr)
            se = min(df.start)-1
            es = max(df.end)+1
        start = Node('start',trs=alltrs,end=se)
        end = Node('end',trs=alltrs,start=es)
        self.start = start
        self.end = end
        self.nodes = {'start':start,'end':end}
        if df is not None:
            for i,row in df.iterrows():
                self.nodes[i] = Node.from_row(i,row)

    def add_node(self,nodename,parents=None, children=None):
        if nodename in self.nodes.keys():
            node = self.nodes[nodename]
        else:
            node = Node(nodename)
            self.nodes[node.name] = node
        if parents is not None:
            for parname in parents:
                if parname not in self.nodes.keys():
                    self.add_node(parname)
                node.add_parent(self.nodes[parname])
                self.nodes[parname].add_child(node)
        if children is not None:
            for childname in children:
                if childname not in self.nodes.keys():
                    self.add_node(childname)
                node.add_child(self.nodes[childname])
                self.nodes[childname].add_parent(node)
    def remove_node(self,name):
        node = self.nodes.pop(name, None)
        if node is not None:
            p = list(node.parents)
            for parent in p:
                parent.remove_child(node)
            c = list(node.children)
            for child in c:
                child.remove_parent(node)
    def __str__(self):
        res = []
        for key,node in self.nodes.items():
            res.append(f"Node {key}, p: [{', '.join(sorted([str(i.name) for i in node.parents]))}], c: [{', '.join(sorted([str(i.name) for i in node.children]))}]")
        return ('\n'.join(res))

    def plot(self):
        dot = Digraph(graph_attr={'rankdir':'LR'})
        for i in self.nodes.keys():
            dot.node(str(i),str(i))
        for key,node1 in self.nodes.items():
            for node2 in node1.children:
                n1,n2=str(key),str(node2.name)
                dot.edge(n1,n2)
        dot
        return dot

def build_graph(df):
    graph = Graph(df = df.copy())
    for i, curr in df.iterrows():
        poss = df.loc[i+1:].copy()
        poss['common'] = poss.tr.apply(lambda x: x&curr.tr)
        poss = poss[poss.common.apply(len)>0]
        children = []
        trs = set(curr.tr)
        for j,child in poss.iterrows():
            if trs&child.tr:
                children.append(j)
                trs-=child.tr
            if not trs:
                break
        else:
            children.append('end')
        graph.add_node(i, children=children)
    for name,node in graph.nodes.items():
        lst = [i.trs for i in node.parents]
        trs = set.union(*lst) if lst else set()
        if (name!='start') and (node.trs-trs):
            node.add_parent(graph.start)
    return graph

def check_pseudoblob(graph,bifurcation):
    for i,row in bifurcation.iterrows():
        for node in graph.nodes.values():
            if (node.name=='start') or (node.name=='end'):
                continue
            if (row.node1.end<node.start) and (node.end<row.node2.start) and \
                                              (set(row.tr)&node.trs) and (node.trs-(set(row.tr))):
                return 1
    return 0

def check_empty_path(graph,bifurcation):
    problem=True
    empty_path = 0
    while problem:
        ind = []
        problem=False
        for i,row in bifurcation.iterrows():
            if row.node2 in row.node1.children: #empty path
                problem=True
                empty_path=1
                ind.append(i)
                #!!!!! we should remove not empty path but all the other paths in this blob
                #row.node1.children = {row.node2}
                #row.node2.parents = {row.node1}
                nodelist = list(graph.nodes.values())
                for node in nodelist:
                    if (node.name=='start') or (node.name=='end'):
                        continue
                    if (node.start>row.node1.end) and (node.end<row.node2.start) and (node.trs&set(row.tr)):
                        graph.remove_node(node.name)
                        
        bifurcation = get_bifurcation(graph)
    return empty_path,bifurcation   

def get_bifurcation(graph):
    starts = [node for node in graph.nodes.values() if len(node.children)>1]
    ends = [node for node in graph.nodes.values() if len(node.parents)>1]
    pairs = []
    for start in starts:
        for end in ends:
            inter = start.trs&end.trs
            if len(inter)>1 and (start.end<end.start):
                pairs.append([start,end,start.end,end.start,tuple(sorted(inter))])
    pairs = pd.DataFrame(pairs,columns=['node1','node2',
                                        'start','end','tr']).sort_values(['start','end']).reset_index(drop=True)
    pairs['l'] = pairs.end-pairs.start
    remove = []
    for i, row1 in pairs.iterrows():
        for j,row2 in pairs.loc[i+1:].iterrows():
            if (row1.end>row2.start) and (row1.tr==row2.tr): #intersects and shares transcript set
                if row1.l>row2.l:
                    remove.append(i)
                else:
                    remove.append(j)
    pairs = pairs[~pairs.index.isin(remove)][['node1','node2','l','tr']]
    lb = [graph.start,graph.end,graph.end.start-graph.start.end,tuple(sorted(graph.start.trs))]
    lb = pd.DataFrame([lb],columns=['node1','node2','l','tr'])
    return pd.concat([pairs,lb]).drop_duplicates(['node1','node2']).sort_values('l')

def reduce_blob(graph,node1,node2,coefs):
    paths = []
    children = [i for i in node1.children if (i.trs&node2.trs) and (i.end<node2.start)]
    trs = set()
    starts = []
    ends = []
    for curr in children:
        path = []
        while curr!=node2:
            path.append(curr.name)
            trs|=curr.trs
            starts.append(curr.start)
            ends.append(curr.end)
            if len(curr.children)>1:
                raise Exception(f'pseudoblob: {curr.name}')
            if len(curr.children)==0:
                raise Exception(f'no children in {curr.name}')
            curr = next(iter(curr.children))
        for name in nmdj.flatten(path):
            coefs[name]/=len(path)
        for name in path:
            graph.remove_node(name)
        paths+=path
    newname = tuple(sorted(nmdj.flatten(paths)))
    newnode = Node(newname,trs=trs,start=min(starts),end=max(ends))
    newnode.parents = {node1}
    newnode.children = {node2}
    node1.children.add(newnode)
    node2.parents.add(newnode)
    graph.nodes[newname] = newnode

def _assign_coefficients(df):
    tags = []
    graph = build_graph(df)
    bifurcation = get_bifurcation(graph)
    pb = check_pseudoblob(graph,bifurcation)
    if pb==1:
        tags.append('pseudo_blob')
        return None,tags
    ep,bifurcation = check_empty_path(graph,bifurcation)
    if ep==1:
        tags.append('empty_path')
    
    coefs = {i:1 for i in graph.nodes.keys() if (i!='start')&(i!='end')}
    for _,row in bifurcation.iterrows():
        reduce_blob(graph,row.node1,row.node2,coefs)
    return coefs,tags



def timeout(seconds_before_timeout):
    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            res = [Exception('function [%s] timeout [%s seconds] exceeded!' % (func.__name__, seconds_before_timeout))]
            def newFunc():
                try:
                    res[0] = func(*args, **kwargs)
                except Exception as e:
                    res[0] = e
            t = Thread(target=newFunc)
            t.daemon = True
            try:
                t.start()
                t.join(seconds_before_timeout)
            except Exception as e:
                print('error starting thread')
                raise e
            ret = res[0]
            if isinstance(ret, BaseException):
                raise ret
            return ret
        return wrapper
    return deco



@timeout(60)
def assign_coefficients(df):
    taglist = []
    try:
        df = df.sort_values(['start','end']).reset_index(drop=True)
        trs = set.union(*df.tr)
        problem = nmdj.check_ir(df)
        ir_flag = False
        while len(problem)>0:
            ir_flag = True
            df = df.drop(problem,axis=0)
            problem = nmdj.check_ir(df)
        trs_new = set.union(*df.tr)
        if trs_new != trs: #deletion of problematic IR leads to inability to quantify some transcripts
            taglist.append('fatal_ir')
            return None,taglist
        nmd = df[df.transcript_biotype=='nonsense_mediated_decay'].copy()
        coding = df[df.transcript_biotype=='protein_coding'].copy()
    
        coding,tags = _assign_coefficients(coding)
        nmd,ntags = _assign_coefficients(nmd)
        tags = ['coding:'+i for i in tags]
        ntags = ['nmd:'+i for i in ntags]
        tags+=ntags
        if ir_flag:
            tags.append('bad_ir')
        taglist+=tags
        if (coding is None) or (nmd is None):
            return None, taglist
        coefs = nmd.copy()
        coefs.update(coding)
    #print(nmd)
    #print(coding)
    #print(coefs)
        df['coef'] = df.index.map(coefs)
        res = df[df.coef.notna()]
        res.loc[res.jtype=='ir','coef'] = res.loc[res.jtype=='ir','coef']/2
        return res, taglist
    except Exception as e:
        print('error in',df.event_id.iloc[0])
        return None, ['error',e]

def get_mat(event):
    trs = set.union(*event.tr)
    mat = np.zeros(shape=(len(trs),event.shape[0]))
    for i,tr in enumerate(trs):
        mat[i:] = event.tr.map(lambda x: tr in x).astype(int).tolist()
    b = np.ones(shape=(len(trs),))
    return mat,b

def nnls_coefficients(event,eps=10e-5):
    try:
        tags = []
        nmd = event[event.transcript_biotype=='nonsense_mediated_decay']
        coding = event[event.transcript_biotype=='protein_coding']
        if nmd.shape[0]>0:
            mat,b = get_mat(nmd)
            ncoef,nnorm = nnls(mat,b)
            if nnorm>eps:
                tags.append('nmd:nnls_norm>0')
                ncoef[:] = np.nan
            nmd['coef'] = ncoef
        else:
            tags.append('nmd:no_junctions')
        if coding.shape[0]>0:
            mat,b = get_mat(coding)
            ccoef,cnorm = nnls(mat,b)
            if cnorm>eps:
                tags.append('coding:nnls_norm>0')
                ccoef[:] = np.nan
            coding['coef'] = ccoef
        else:
            tags.append('coding:no_junctions')
        res = pd.concat([nmd,coding])
        res.loc[res.jtype=='ir','coef'] = res.loc[res.jtype=='ir','coef']/2
        res = res[res.coef.notna()]
        if res.transcript_biotype.nunique()==2:
            return res,tags
        return None,tags
    except Exception as e:
        print('error in',event.event_id.iloc[0])
        return None, ['error',e]

def get_coefs(results,threads=5,start_time=None,use_nnls=False):
    tag_list = []
    output = []
    if start_time is None:
        start_time = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        parallel_func = nnls_coefficients if use_nnls else assign_coefficients
        futures = dict()
        for event, temp in results.groupby('event_id'):
            futures[executor.submit(parallel_func, temp)] = event
        print(f'tasks created, time: {round((time.time()-start_time)/60, 2)} min')
        for i,future in enumerate(concurrent.futures.as_completed(futures,timeout=3*60*60)):
                try:
                    coef,tags = future.result()
                except Exception as e:
                    print(futures[future])
                    print(e)
                    tag_list.append([futures[future],('parallel_error',e)])
                    continue
                tag_list.append([futures[future],tags])
                if coef is not None:
                    output.append(coef)
                if i%1000==0:
                    print(f'Done {i} events, time: {time.time()-start_time}')
    if output:
        output = pd.concat(output)
    else:
        output=None
    tags = pd.DataFrame(tag_list,columns=['event_id','coef_tags'])
    #tagu = set(nmdj.flatten(tags.coef_tags))
    #for t in tagu:
    #    tags[str(t)] = tags.coef_tags.apply(lambda x: t in x).astype('int')
    return output,tags