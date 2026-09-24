"""Resume the immutable screen with retries for transient Windows/OneDrive locks."""
import os,time,json
from pathlib import Path
import tone_count_study as screen

def resilient_save(path,x):
 path.parent.mkdir(parents=True,exist_ok=True)
 tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
 tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
 for attempt in range(100):
  try:tmp.replace(path);return
  except PermissionError:
   if attempt==99:raise
   time.sleep(.1)

if __name__=='__main__':
 screen.save=resilient_save
 screen.main()
