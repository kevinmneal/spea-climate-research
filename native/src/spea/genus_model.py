"""Small native BIO19 reader and stratified quadrature for genus development."""
from pathlib import Path
import numpy as np
from affine import Affine
from pyproj import Transformer
from .genus_spatial import ALBERS

LOG_KEYS={'bio12','bio13','bio14','bio16','bio17','bio18','bio19'}

class NativeBIO19:
    def __init__(self,root,registry):
        self.fields=registry['fields'];self.keys=[f'bio{i:02}' for i in range(1,20)]
        by_key={x['key']:x for x in self.fields}
        self.arrays=[np.load(Path(root)/by_key[k]['physical_path'],mmap_mode='r').reshape(-1) for k in self.keys]
        self.shape=tuple(registry['window_grid']['shape'])
        if any(len(x)!=np.prod(self.shape) for x in self.arrays):raise ValueError('Native feature grid mismatch')

    def read(self,indices):
        ids=np.asarray(indices,dtype=np.int64)
        if ids.ndim!=1 or np.any(ids<0) or np.any(ids>=np.prod(self.shape)):raise ValueError('Invalid native flat indices')
        x=np.column_stack([a[ids] for a in self.arrays]).astype('float64')
        if not np.isfinite(x).all():raise ValueError('Missing native climate; no filling or snapping')
        for j,k in enumerate(self.keys):
            if k in LOG_KEYS:
                if np.any(x[:,j]<0):raise ValueError('Negative precipitation amount')
                x[:,j]=np.log1p(x[:,j])
        return x

def geographic_blocks(flat_indices,shape,transform,width_m):
    ids=np.asarray(flat_indices,dtype=np.int64);row,col=np.divmod(ids,shape[1])
    lon,lat=Affine(*transform[:6])@(col+.5,row+.5)
    x,y=Transformer.from_crs(4326,ALBERS,always_xy=True).transform(lon,lat)
    bx=np.floor(np.asarray(x)/width_m).astype('int64');by=np.floor(np.asarray(y)/width_m).astype('int64')
    # Injective within explicitly checked regional grid; no biological meaning.
    if np.any(np.abs(by)>=5000):raise ValueError('Block encoding outside supported regional domain')
    return bx*10000+by

def stratified_quadrature(eligible,forced,block_ids,areas,per_block,seed):
    """Stable sorted O(N log N) equivalent of native_benchmark quadrature_sample.

    block_ids and areas align with eligible, rather than an entire global grid.
    Return native IDs, integration weights, inclusion probabilities and strata.
    """
    eligible=np.asarray(eligible,dtype=np.int64);forced=np.asarray(forced,dtype=np.int64)
    block_ids=np.asarray(block_ids);areas=np.asarray(areas,dtype=float)
    if len(np.unique(eligible))!=len(eligible) or len(np.unique(forced))!=len(forced):raise ValueError('Unique quadrature and presence cells required')
    if not np.isin(forced,eligible).all():raise ValueError('Every presence must be in the quadrature domain')
    if block_ids.shape!=eligible.shape or areas.shape!=eligible.shape or per_block<1 or not np.isfinite(areas).all() or np.any(areas<=0):raise ValueError('Invalid quadrature strata or areas')
    isforced=np.isin(eligible,forced);rest=np.flatnonzero(~isforced)
    order=np.argsort(block_ids[rest],kind='stable');ordered=rest[order]
    cuts=np.r_[0,np.flatnonzero(np.diff(block_ids[ordered]))+1,len(ordered)]
    selected=[np.flatnonzero(isforced)];probability=[np.ones(isforced.sum())];counts=[]
    rng=np.random.default_rng(seed)
    for a,b in zip(cuts[:-1],cuts[1:]):
        pool=ordered[a:b]
        if not len(pool):continue
        size=min(per_block,len(pool));take=rng.choice(pool,size,replace=False)
        selected.append(take);probability.append(np.full(size,size/len(pool)))
        counts.append((int(block_ids[pool[0]]),len(pool),size))
    selected=np.concatenate(selected);pi=np.concatenate(probability);order=np.argsort(eligible[selected]);selected=selected[order];pi=pi[order]
    return eligible[selected],areas[selected]/pi,pi,block_ids[selected],np.asarray(counts,dtype=np.int64).reshape(-1,3)
