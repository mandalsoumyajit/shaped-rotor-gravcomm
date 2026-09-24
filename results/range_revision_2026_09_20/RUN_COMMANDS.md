# Calculation commands and provenance

Executed on Windows using C:/ProgramData/anaconda3/python.exe (Python 3.12;
NumPy 1.26.4, SciPy 1.13.1, Numba 0.60.0). Scripts set BLAS/OpenMP threads to one.
The commands below are relative to AIP_Advances_Submission; `python` denotes that runtime.

```text
python scripts/range_revision.py prepare
python scripts/range_revision.py screen --packets 32 --workers 4
python scripts/range_revision.py validate --cases desktop_2m --max-packets 10000 --workers 4
python scripts/range_revision.py validate --max-packets 25000 --workers 4
python scripts/range_array_comparison.py
python scripts/range_revision_checks.py --packets 1000 --workers 4
python tests/test_range_revision.py
python scripts/range_revision.py validate --cases desktop_1.5m --symbol-s 1259.93 --q 0.5 --tag followup_q05 --alpha 0.025 --max-packets 15000 --workers 4
python scripts/range_revision.py validate --cases scale2_2m_same_receiver --symbol-s 25.3833 --q 0.5 --tag conservative --alpha 0.025 --max-packets 10000 --workers 4
python scripts/range_revision.py validate --cases scale2_2m_slow_retuned_receiver --symbol-s 14.3602 --q 0.5 --tag conservative --alpha 0.025 --max-packets 10000 --workers 4
python scripts/range_revision.py validate --cases desktop_2m --max-packets 40000 --workers 2
python scripts/range_revision.py validate --cases desktop_1.5m --symbol-s 1259.93 --q 0.5 --tag followup_q05 --alpha 0.025 --max-packets 25000 --workers 2
python scripts/range_revision_checks.py --cases scale2_2m_same_receiver scale2_2m_slow_retuned_receiver --tag conservative --from-validation --packets 1000 --workers 4
python scripts/summarize_range_revision.py
```

Validation streams resume by packet index. Different output tags use independent
random seeds. Do not change candidate parameters or confidence settings within
a resumed file; the script checks them. Source quadrature SHA-256 is recorded in
scenarios.json. Exact candidates, error counts per packet, and seeds are in JSON.

The original coarse-grid candidate is retained even when it failed or remained
inconclusive. The report recomputes the packet-confidence interval at alpha=0.05/K,
where K is the number of independently validated candidate files for that scenario.
This protects the within-scenario choice; it is not a simultaneous statement across
all scenarios. New future candidate trials require updating K and rechecking bounds.

Status: initial calculations, not final optimization. New actuator trajectories,
environmental-noise sensitivity, physical acquisition/phase drift, and an independent
parallel-subchannel array simulation remain to be completed. The accepted manuscript,
baseline BER data, and earlier scaled-source results have not been overwritten.
