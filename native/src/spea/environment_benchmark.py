"""Matched presence/background feature comparisons using ordinary classifiers.

Raw classifier log odds estimate an occurrence/background density ratio under
the declared weighted case-control sample. They are not occupancy probabilities.
No conditional point-process gradient is passed through a classifier API.
"""
from pathlib import Path
import numpy as np
import pandas as pd


def case_control_weights(n_presence, background_area):
    """Equal total weight for each class; physical area within background.

    Mean object weight is one. This scale is kept identical between predictor
    families for a given split; an arbitrary change in weight scale can change
    a regularized fit. Background area units cancel.
    """
    a=np.asarray(background_area,dtype=float)
    if n_presence<1 or a.ndim!=1 or not len(a) or not np.isfinite(a).all() or np.any(a<=0):
        raise ValueError('Require positive presence count and finite positive background weights')
    half=(n_presence+len(a))/2
    scaled=a/a.max();scaled/=scaled.sum()
    return np.r_[np.full(n_presence,half/n_presence),half*scaled]


class NativeEnvironment:
    """Chunk reader that closes each memory map before opening the next one."""
    def __init__(self,root,shape,fields):
        self.root=Path(root);self.shape=tuple(shape)
        self.fields={f['key']:f for f in fields}
        if len(self.fields)!=len(fields):raise ValueError('Repeated environmental key')

    def read(self,indices,keys):
        ids=np.asarray(indices,dtype=np.int64)
        if ids.ndim!=1 or np.any(ids<0) or np.any(ids>=np.prod(self.shape)):
            raise ValueError('Indices outside native window')
        data={}
        for key in keys:
            field=self.fields[key]
            if 'missing_of' in field:
                value=self.read(ids,field['missing_of']).isna().any(axis=1).astype('float32').to_numpy()
            else:
                array=np.load(self.root/field['physical_path'],mmap_mode='r')
                if tuple(array.shape)!=self.shape:raise ValueError('Environmental grid shape mismatch: '+key)
                value=np.array(array.reshape(-1)[ids],dtype='float32',copy=True)
                array._mmap.close()
                if not field.get('allow_missing',False) and not np.isfinite(value).all():
                    raise ValueError('Unexpected missing environmental values: '+key)
                if np.isinf(value).any():raise ValueError('Infinite environmental value: '+key)
                if field.get('transform')=='log1p':
                    if np.any(value[np.isfinite(value)]<0):raise ValueError('Negative amount: '+key)
                    value=np.log1p(value)
                if field.get('categorical',False):
                    if not np.isfinite(value).all() or np.any(value!=np.floor(value)):
                        raise ValueError('Categorical IDs must be finite integers')
                    value=value.astype('int64').astype(str)
            data[key]=value
        return pd.DataFrame(data)
