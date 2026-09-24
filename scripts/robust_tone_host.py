"""I/O-resilient host wrapper; simulation code/manifests remain unchanged."""
import sys,importlib
from resume_tone_count_screen import resilient_save
if __name__=='__main__':
 name={'gates':'tone_count_gates','confirm':'tone_count_confirm','extra':'tone_count_extra'}[sys.argv[1]]
 module=importlib.import_module(name)
 if name=='tone_count_gates':module.runtime.save=resilient_save
 else:module.save=resilient_save
 module.main()
