#!/bin/bash
# Final hop-limited CPFSK FEA run, tag 2026_11_05_nb. Run in WSL Ubuntu-24.04 from a Linux-side copy of the release
# (the solver copies every output back into fea/results/*_2026_11_05_nb).
# Environment: CalculiX 2.21 (Ubuntu calculix-ccx 2.21-1, user-local), Ubuntu python3-gmsh 4.12.1+ds1-1.1build2
# (user-local) and Ubuntu system python3/numpy/scipy, i.e. the archived 2026-09-20 package set. The PyPI gmsh 4.12.1
# wheel does NOT reproduce the archived tungsten meshes byte-for-byte; the Ubuntu build does.
set -euo pipefail
source $HOME/local/gravcomm_fea_env_deb.sh
PY=/usr/bin/python3
export GRAVCOMM_FEA_TAG=2026_11_05_nb
cd "${FEA_DIR:-$HOME/gravcomm_fea_nb_2026_11_05/public_release_v0_5_0/fea}"
S=results/three_link_designs_$GRAVCOMM_FEA_TAG; T=results/tungsten_comparison_$GRAVCOMM_FEA_TAG
stamp(){ echo "=== $(date -Is) $*"; }
stamp versions; { ccx -v 2>&1 || true; } | grep -i version; $PY -c "import gmsh,numpy,scipy,sys;print(sys.version.split()[0],gmsh.__version__,gmsh.__file__,numpy.__version__,scipy.__version__)"
stamp frozen designs;         $PY -s make_hop_frozen_designs.py > /dev/null
stamp steel solve;            $PY -s three_link_fea_hop.py --prepare --workers 8
stamp steel compile .018;     $PY -s compile_three_link_fea_hop.py --mesh .018
stamp steel compile .013;     $PY -s compile_three_link_fea_hop.py --mesh .013
stamp steel modal .018;       $PY -s three_link_modal_response_hop.py --mesh .018
stamp steel modal .013;       $PY -s three_link_modal_response_hop.py --mesh .013
stamp assembly;               $PY -s assembly_screen_hop.py > $S/assembly_screen.log
stamp steel report;           $PY -s report_three_link_fea_hop.py
stamp tungsten solve;         $PY -s tungsten_link_fea_hop.py --prepare --workers 8
stamp tungsten compile .018;  $PY -s compile_tungsten_link_fea_hop.py --mesh .018
stamp tungsten compile .013;  $PY -s compile_tungsten_link_fea_hop.py --mesh .013
stamp tungsten modal .018;    $PY -s tungsten_modal_response_hop.py --mesh .018
stamp tungsten modal .013;    $PY -s tungsten_modal_response_hop.py --mesh .013
stamp tungsten field;         $PY -s tungsten_field_comparison_hop.py
stamp tungsten refine h0.010; $PY -s refine_tungsten_reference_hop.py
stamp tungsten report;        $PY -s report_tungsten_comparison_hop.py
stamp summary;                $PY -s summarize_hop_vs_archived.py > $S/hop_vs_archived_summary.log
stamp done
