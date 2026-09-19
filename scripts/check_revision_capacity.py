"""Band/grid convergence using the shaped rotor's stored finite-volume field."""
import csv,json
from pathlib import Path
import numpy as np
from paper_results import gaussian_waterfill
from gravcomm.receivers import CARTER_OSCILLATOR as receiver

ROOT=Path(__file__).resolve().parents[1]
sources=[r for r in csv.DictReader((ROOT/'fea/results/shaped_waveforms/summary.csv').open())
         if r['case']=='medium_rounded']
records=[]
for source in sources:
    amplitude=float(source['a2_m_s2'])
    for bandwidth in [1.5625,2.,4.,8.]:
        values=[]
        for n in [32768,65536]:
            f=50.3+bandwidth*((np.arange(n)+.5)/n-.5)
            values.append(gaussian_waterfill(receiver.acceleration_noise_asd(f)**2,
                                            bandwidth,amplitude**2/2))
        records.append(dict(distance_m=float(source['distance_m']),amplitude_m_s2=amplitude,
            bandwidth_hz=bandwidth,capacity_bit_s=values[1],
            relative_grid_change=values[1]/values[0]-1))
assert max(abs(r['relative_grid_change']) for r in records)<1e-5
(ROOT/'results/communications_checks/capacity.json').write_text(json.dumps(records,indent=2)+'\n')
noise_records=[]
amplitude=float(next(r for r in sources if float(r['distance_m'])==.5)['a2_m_s2'])
for added in [0.,1e-10,2e-10]:
    for bandwidth in [2.,4.,8.]:
        f=50.3+bandwidth*((np.arange(65536)+.5)/65536-.5)
        noise_records.append(dict(added_acceleration_asd=added,bandwidth_hz=bandwidth,
            capacity_bit_s=gaussian_waterfill(receiver.acceleration_noise_asd(f)**2+added**2,
                                              bandwidth,amplitude**2/2)))
(ROOT/'results/communications_checks/technical_noise_capacity.json').write_text(json.dumps(noise_records,indent=2)+'\n')
print(json.dumps(records,indent=2))
