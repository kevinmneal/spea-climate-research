"""Fixed geographic fold assignments and conservative WGS84 separation bounds.

No environmental values or occurrence outcomes enter the lattice. Distances use
ellipsoid bounds or exact geodesics, never Albers map units as exact distances.
"""
from __future__ import annotations
import numpy as np
from affine import Affine
from pyproj import Geod,Transformer
from scipy.ndimage import binary_erosion
from scipy.spatial import cKDTree

ALBERS='+proj=aea +lat_1=20 +lat_2=50 +lat_0=35 +lon_0=-105 +datum=WGS84 +units=m +no_defs'
GEOD=Geod(ellps='WGS84')
A=6378137.;F=1/298.257223563;E2=F*(2-F)
MIN_METRIC_RADIUS=A*(1-E2)
MAX_METRIC_RADIUS=A/np.sqrt(1-E2)

def sphere_xyz(longitude,latitude):
 lon,lat=np.deg2rad(longitude),np.deg2rad(latitude)
 if not np.isfinite(lon+lat).all() or np.any(np.abs(lat)>np.pi/2):raise ValueError('Finite valid coordinates required')
 return np.column_stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])

def lattice(longitude,latitude,width_m=150000,folds=5):
 if width_m<=0 or folds!=5:raise ValueError('This frozen five-colour design requires positive width and five folds')
 x,y=Transformer.from_crs(4326,ALBERS,always_xy=True).transform(longitude,latitude)
 ix=np.floor(np.asarray(x)/width_m).astype('int32');iy=np.floor(np.asarray(y)/width_m).astype('int32')
 return ix,iy,((ix+2*iy)%folds).astype('int8')

def pixel_cover_radius(transform):
 """Global upper bound from center to any point of an angular raster cell."""
 t=Affine(*transform[:6])
 if t.b!=0 or t.d!=0 or t.a<=0 or t.e>=0:raise ValueError('North-up angular raster required')
 return MAX_METRIC_RADIUS*np.deg2rad((abs(t.a)+abs(t.e))/2)

def rectangle_inside_clearance(longitude,latitude,bounds_wsen):
 """Conservative WGS84 clearance to an angular rectangle's exterior.

 Latitude change bounds distance to each parallel. Distance to the complete
 meridian great circle bounds distance to its finite boundary segment. Scaling
 unit-sphere distances by the global minimum ellipsoid metric radius preserves
 both lower bounds; taking their minimum covers all four edges. This avoids a
 slow nearest-neighbour search along an extremely thin rectangular frontier.
 """
 west,south,east,north=bounds_wsen
 if not (-180<=west<east<=180 and east-west<180 and -90<south<north<90):raise ValueError('Non-wrapping rectangle under180degrees required')
 lon,lat=np.asarray(longitude),np.asarray(latitude);phi=np.deg2rad(lat)
 ns=np.minimum(np.abs(phi-np.deg2rad(south)),np.abs(phi-np.deg2rad(north)))
 ew=np.minimum(np.arcsin(np.clip(np.abs(np.cos(phi)*np.sin(np.deg2rad(lon-west))),0,1)),
               np.arcsin(np.clip(np.abs(np.cos(phi)*np.sin(np.deg2rad(east-lon))),0,1)))
 result=MIN_METRIC_RADIUS*np.minimum(ns,ew)
 result[(lon<west)|(lon>east)|(lat<south)|(lat>north)]=-np.inf
 return result

class RasterFrontier:
 """Distance lower bounds to the edge of a held geographic pixel footprint.

 All true-mask border cells, including hole/cache borders, cover the boundary.
 WGS84 surface distance >= MIN_METRIC_RADIUS * unit-sphere angular distance.
 Subtracting the full angular-cell cover bound gives a conservative boundary
 clearance. Direct mask membership is checked separately: an interior point
 cannot be mistaken for safe training merely because it is far from the edge.
 """
 def __init__(self,mask,transform):
  self.mask=np.asarray(mask,dtype=bool);self.transform=Affine(*transform[:6]);self.cover=pixel_cover_radius(transform)
  if self.mask.ndim!=2 or not self.mask.any():raise ValueError('Nonempty two-dimensional domain required')
  border=self.mask&~binary_erosion(self.mask,structure=np.ones((3,3),bool),border_value=0)
  self.rows,self.cols=np.nonzero(border)
  self.longitude,self.latitude=self.transform@(self.cols+.5,self.rows+.5)
  self.tree=cKDTree(sphere_xyz(self.longitude,self.latitude))
 def contains(self,longitude,latitude):
  col,row=(~self.transform)@(np.asarray(longitude),np.asarray(latitude))
  row=np.floor(row).astype('int64');col=np.floor(col).astype('int64')
  inside=(row>=0)&(col>=0)&(row<self.mask.shape[0])&(col<self.mask.shape[1])
  result=np.zeros(inside.shape,bool);result[inside]=self.mask[row[inside],col[inside]]
  return result
 def boundary_lower_bound(self,longitude,latitude,chunk=65536):
  lon,lat=np.asarray(longitude),np.asarray(latitude);result=np.empty(len(lon),float)
  for start in range(0,len(lon),chunk):
   stop=min(start+chunk,len(lon))
   chord,_=self.tree.query(sphere_xyz(lon[start:stop],lat[start:stop]),workers=1)
   angle=2*np.arcsin(np.clip(chord/2,0,1))
   result[start:stop]=MIN_METRIC_RADIUS*angle-self.cover
  return result
 def outside_clearance(self,longitude,latitude,query_radius_m=0):
  result=self.boundary_lower_bound(longitude,latitude)-np.asarray(query_radius_m)
  result[self.contains(longitude,latitude)]=-np.inf
  return result
 def inside_clearance(self,longitude,latitude,support_radius_m=0):
  result=self.boundary_lower_bound(longitude,latitude)-np.asarray(support_radius_m)
  result[~self.contains(longitude,latitude)]=-np.inf
  return result

def cell_support_radius(records,cell_longitude,cell_latitude,unknown_radius_m):
 """Enclose every contributor's explicit/assumed disk about the cell center.

 The offset is geometry, not added provider precision. Missing source radii
 remain unchanged; the returned value belongs only to a named scenario.
 """
 if unknown_radius_m<0:raise ValueError('Nonnegative scenario radius required')
 known=records.coordinate_uncertainty_m.to_numpy(float)
 if np.isinf(known).any() or np.any(known<0):raise ValueError('Invalid reported support radius; do not replace with unknown scenario')
 chosen=np.where(np.isfinite(known)&(known>0),known,unknown_radius_m)
 offsets=np.abs(GEOD.inv(np.full(len(records),cell_longitude),np.full(len(records),cell_latitude),records.longitude.to_numpy(),records.latitude.to_numpy())[2])
 return float(np.max(chosen+offsets))

def exact_pair_min_clearance(train_lon,train_lat,train_radius,test_lon,test_lat,test_radius):
 """Exact WGS84 center distance minus both enclosing support radii."""
 tlon,tlat,tr=np.asarray(train_lon),np.asarray(train_lat),np.asarray(train_radius)
 result=np.full(len(tlon),np.inf)
 for lon,lat,radius in zip(test_lon,test_lat,test_radius):
  d=np.abs(GEOD.inv(tlon,tlat,np.full(len(tlon),lon),np.full(len(tlat),lat))[2])
  result=np.minimum(result,d-tr-radius)
 return result
