"""Historical nominal Macnae magnetic-background ASD, manually digitized.

Returns T/sqrt(Hz); simulation PSD is ASD squared (one-sided convention).
No extrapolation or interpolation across the below-axis gap is permitted.
This schematic does not specify line bandwidths, site statistics or covariance.
"""
from pathlib import Path
from functools import lru_cache
import csv
import numpy as np
DEFAULT_DATA=Path(__file__).resolve().parents[2]/'results/macnae_noise_digitization_2026_09_21/anchors.csv'

@lru_cache(None)
def _anchors(path):
 with Path(path).open() as stream:rows=list(csv.DictReader(stream))
 return tuple((np.array([float(r['frequency_Hz']) for r in rows if r['segment']==seg]),np.array([float(r['asd_T_sqrtHz']) for r in rows if r['segment']==seg])) for seg in ('low','high'))

def macnae_asd(frequency_Hz,scale=1.,unsupported='raise',data_path=DEFAULT_DATA):
 """Piecewise power-law interpolation; scale multiplies field ASD.

 unsupported='nan' marks untraced frequencies; default raises ValueError.
 Narrow spectral lines require separate user-specified RMS and linewidth.
 """
 f=np.asarray(frequency_Hz,dtype=float)
 if np.any(~np.isfinite(f)) or np.any(f<=0):raise ValueError('Frequencies must be finite and positive')
 if not np.isfinite(scale) or scale<=0:raise ValueError('ASD scale must be finite and positive')
 if unsupported not in ('raise','nan'):raise ValueError('unsupported must be raise or nan')
 out=np.full(f.shape,np.nan)
 for x,y in _anchors(str(data_path)):
  mask=(f>=x[0])&(f<=x[-1])
  out=np.where(mask,np.exp(np.interp(np.log(f),np.log(x),np.log(y)))*scale,out)
 if unsupported=='raise' and np.any(np.isnan(out)):raise ValueError('Frequency outside traced segments or in unresolved gap; use unsupported="nan" to mark it')
 return out
