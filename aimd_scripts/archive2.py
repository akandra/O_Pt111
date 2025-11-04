#%%
from pathlib import Path
import tarfile

for d in Path().glob('*'):
    if d.is_dir():
        try:
            int(d.name)
            if int(d.name)<3:
                with tarfile.open(d.name+".tar.gz", "w:gz") as tar: tar.add(d)
        except: pass




