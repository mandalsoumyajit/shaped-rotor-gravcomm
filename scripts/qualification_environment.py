import sys,json,platform,os
import numpy,scipy,numba
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p=root/'results/fixed_rate_sizing_2026_09_20/ber_qualification_v1/environment.json'
p.write_text(json.dumps(dict(python=sys.version,executable=sys.executable,platform=platform.platform(),machine=platform.machine(),numpy=numpy.__version__,scipy=scipy.__version__,numba=numba.__version__,logical_cpus=os.cpu_count(),worker_count=40,blas_thread_settings='streaming_dynamics_runtime sets OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, OMP_NUM_THREADS to 1 before importing numpy'),indent=2)+'\n')
print(p.read_text())
