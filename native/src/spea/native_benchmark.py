"""Small native-grid benchmark helpers, shared by preparation and fitting."""
import json
from pathlib import Path
import numpy as np
from pyproj import Geod
from scipy.special import logsumexp

GEOD=Geod(ellps='WGS84')


def features_for_arm(directory, arm, indices):
    directory=Path(directory)
    key='bio19' if arm=='bio19' else 'monthly27'
    result=np.load(directory/f'{key}.npy',mmap_mode='r')[indices]
    if arm=='monthly29_soil':
        result=np.column_stack([result,np.load(directory/'soil2.npy',mmap_mode='r')[indices]])
    if not np.isfinite(result).all():
        raise ValueError(f'Incomplete features for {arm}')
    return result


def region_distance_m(latitude, breaks, region):
    """Meridional distance to the nearest boundary of a latitude-band region."""
    latitude=np.asarray(latitude)
    lo=[-90,*breaks][region];hi=[*breaks,90][region]
    edge=np.clip(latitude,lo,hi)
    return np.abs(GEOD.inv(np.zeros(len(latitude)),latitude,np.zeros(len(latitude)),edge)[2])


def support_guard(train, held, longitude, latitude, radius, buffer_m):
    """Remove training supports too close to any held support; never alter data."""
    train=np.asarray(train,dtype=int);held=np.asarray(held,dtype=int)
    if not len(train) or not len(held):
        raise ValueError('Nonempty training and held presence sets required')
    t=np.repeat(train,len(held));h=np.tile(held,len(train))
    distances=np.abs(GEOD.inv(longitude[t],latitude[t],longitude[h],latitude[h])[2]).reshape(len(train),len(held))
    separation=distances-radius[train,None]-radius[held][None,:]
    keep=separation.min(axis=1)>=buffer_m
    return train[keep],train[~keep],float(separation[keep].min()) if keep.any() else None


def quadrature_sample(eligible, forced, blocks, area, per_block, seed):
    """Stratified without-replacement quadrature with explicit inclusion weights."""
    eligible=np.asarray(eligible,dtype=int);forced=np.asarray(forced,dtype=int)
    if not np.isin(forced,eligible).all() or len(np.unique(forced))!=len(forced):
        raise ValueError('Forced presence cells must be unique and inside quadrature domain')
    if per_block<1 or not np.isfinite(area[eligible]).all() or np.any(area[eligible]<=0):
        raise ValueError('Positive areas and at least one sample per stratum required')
    rng=np.random.default_rng(seed)
    rest=eligible[~np.isin(eligible,forced)]
    selected=[forced];probabilities=[np.ones(len(forced))]
    for block in np.unique(blocks[rest]):
        pool=rest[blocks[rest]==block];count=min(per_block,len(pool))
        selected.append(rng.choice(pool,count,replace=False))
        probabilities.append(np.full(count,count/len(pool)))
    selected=np.concatenate(selected);pi=np.concatenate(probabilities)
    order=np.argsort(selected)
    return selected[order],(area[selected]/pi)[order],pi[order]


def predict_chunked(model, directory, arm, indices, chunk_size=16384):
    values=[]
    for start in range(0,len(indices),chunk_size):
        values.append(model.predict_log_intensity(features_for_arm(directory,arm,indices[start:start+chunk_size])))
    return np.concatenate(values) if values else np.empty(0)


def relative_scores(presence_eta, background_eta, area):
    if not len(background_eta) or not np.isfinite(background_eta).all():
        raise ValueError('Finite target integration required')
    normalizer=logsumexp(background_eta+np.log(area))-np.log(area.sum())
    return np.asarray(presence_eta)-normalizer,float(normalizer)


def block_score(values, blocks):
    unique,inverse=np.unique(blocks,return_inverse=True)
    if not len(unique):return float('nan')
    return float(np.mean(np.bincount(inverse,weights=values)/np.bincount(inverse)))
