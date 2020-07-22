import numpy as np 
import shutil
import sys
import os 
import glob 

NEK5K  = sys.argv[1]  
INSITU = {}  
INSITU["V"] = (NEK5K,5)
INSITU["W"] = (NEK5K,6)
INSITU["X"] = (NEK5K,9)  
INSITU["Y"] = (NEK5K,8)   
INSITU["Z"] = (NEK5K,7)   

##
for k,v in INSITU.items() :
  files = glob.glob('./%s*' % k )
  for f in files : os.system("rm -rf %s" % f)
  os.symlink("pipe.py", "%s.py"%k)


## 'mpmd.conf'
CumSum = np.cumsum([0]+[v[1] for v in INSITU.values()])
Bottom = dict(zip(INSITU.keys(),CumSum[ :-1]))
Top    = dict(zip(INSITU.keys(),CumSum[1:  ]))
SMS    = "\n".join([" %d-%d %s --namei %s " % (Bottom[k],Top[k]-1,v[0],k) for k,v in INSITU.items()]) 

#print( SMS )
open("mpmd.conf","w").write(SMS)

