"""Run the paper-study checks (paper_studies/ground_motion_study/test_joint_receiver.py)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'paper_studies'/'ground_motion_study'))
from test_joint_receiver import *  # noqa: F401,F403
