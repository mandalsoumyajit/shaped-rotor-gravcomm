# Release validation

Validation performed in the independent v0.4.0 release checkout on 2026-09-21:

- 117 numerical and communications tests passed (Python/NumPy/SciPy/Numba).
- Five FEA parser, stress, shape-function and averaging tests passed in WSL with Gmsh installed.
- All scientific source hashes match both frozen BER qualification manifests.
- The release data archive retains the original packet outcomes and numerical inputs; serialized workstation folder paths are repository-relative.

The release adds no new Monte Carlo samples or hardware claims. The numerical conclusions and confidence-family limitations are documented in the included results reports. Packet audits and final artifact integrity checks are recorded below after completion.

- Both full qualification compilers passed: all 36,992 records, source hashes, detector gates and first stopping prefixes verified.
- Steel and tungsten mechanical source hashes match their respective frozen manifests.
- Current summary figure regenerated successfully from the release checkout.
- All 39,081 data-archive entries verified against their SHA-256 manifest after path relocation.
