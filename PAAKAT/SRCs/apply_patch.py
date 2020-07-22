import os 
import sys 
import argparse
os.system('clear')

#import md5
def CheckSum( fname ):
  hash1 = md5.new()
  hash1.update( fname )
  return hash1.digest()

import hashlib
def CheckSum( fname ):
  hasher = hashlib.md5()
  with open(fname, 'rb') as afile:
    buf = afile.read()
    hasher.update(buf)
  #print(hasher.hexdigest())
  return hasher.hexdigest() 


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-O', action='store', dest='Origin',
                    default=[], type=str, nargs='+',
                    help='Orig_Path',
                   )
parser.add_argument('-D', action='store', dest='Destination',
                    default=[], type=str, nargs='+',
                    help='Destination_Path',
                   )
options = parser.parse_args()
#
try :
  options.Origin[0] 
  options.Destination[0]  
except : 
  parser.print_help()
  print() 
  sys.exit()

options = vars(options)
for k,v in options.items() : print("\t[apply_patch] k:'%s', v:" % k,v)


# -------------------------------------------------------------------------
GetFiles    = lambda _path : [(R,f) for R,D,F in os.walk(_path) for f in F]
GetFilesDic = lambda _path : {f:R for R,D,F in os.walk(_path) for f in F} 


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
OriginPath  = options['Origin'][0]
OriginFiles = GetFilesDic(OriginPath) ; 
#for f in OriginFiles : print(f) 
DestPath    = options['Destination'][0]
DestFiles   = GetFiles(DestPath) ; ##for f in DestFiles : print(f)  

import shutil
import filecmp
from datetime import datetime
now = datetime.now().strftime('%Y%B%d_%H%M%S')

Founds = []  
#for src_name,src_path in OriginFiles.items() : 
DictChangedSizeDuringIter = {k:v for k,v in OriginFiles.items() if v}  
for src_name,src_path in DictChangedSizeDuringIter.items() :
  src_orig = os.path.join(src_path,src_name); #print("\t[apply_patch] SrcOrig:'%s'" % src_orig )
  os.path.getmtime( src_orig );   
  CheckSum(src_orig) 
  for D in DestFiles :   
    if(D[1]==src_name) :
      #print( D[0] )  
      #Founds.append( D )

      if 1 : 
        if OriginFiles.get(src_name,0) : del OriginFiles[src_name]  

        ## 'symlink' style. SEE [1] 
        dest_name = src_orig.replace(OriginPath,''); #print( dest_name )
        dest_orig = DestPath + dest_name; #print("\t[apply_patch] DstOrig:'%s'" % dest_orig )
        dest_path,dest_name = os.path.split(dest_orig); #print( dest_path,dest_name,src_name ) 
        assert( os.path.exists(dest_path) )

        if os.path.exists(dest_orig) :
          #print("\t '%s' exist!" % dest_orig) 
          dest_copy = os.path.join(dest_path,"%s_%s"%(dest_name,now));
          equal     = filecmp.cmp(src_orig, dest_orig);
          if equal :
            print("\t[apply_patch] '%s' NO CHANGES ('%s') " %(dest_name,D[0]) )
          else :
            print("\t[apply_patch] '%s' -> '%s' PATCHED " %(dest_name,"%s_%s"%(dest_name,now)) )
            shutil.move(dest_orig, dest_copy)
            os.symlink(src_orig, dest_orig) 
          Founds.append( D )
        else :
          #print("\t[apply_patch] '%s' no found!" % dest_orig)
          os.symlink(src_orig, dest_orig)
          print("\t[apply_patch] '%s' NEW ('%s') " %(src_name,D[0]) )
          Founds.append( D )

      else :
        ## 'original' style  (first implementation) 
        dest_path = D[0] # path 
        dest_name = D[1] # name 
        dest_orig = os.path.join(dest_path,dest_name); #print("\t[apply_patch] DstOrig:'%s'" % dest_orig )
      
        dest_copy = os.path.join(dest_path,"%s_%s"%(dest_name,now)); #print("\t[apply_patch] DstCopy:'%s'" % dest_copy )
        equal     = filecmp.cmp(src_orig, dest_orig); ##print("\t[apply_patch] %s " % equal) 
        if equal : 
          print("\t[apply_patch] '%s' no changes " %(dest_name) )
        else : 
          print("\t[apply_patch] '%s' -> '%s' patched " %(dest_name,"%s_%s"%(dest_name,now)) )   
          shutil.move(dest_orig, dest_copy)
          os.symlink(src_orig, dest_orig)  

#exit(0)
## Not found (New) files in 'DestPath'
#for src_name,src_path in OriginFiles.items() : 
DictChangedSizeDuringIter = {k:v for k,v in OriginFiles.items() if v}
for src_name,src_path in DictChangedSizeDuringIter.items() :
  src_orig  = os.path.join(src_path,src_name); print("\t[apply_patch] SrcOrig:'%s'" % src_orig )

  ## [1] 'symlink' style. 
  dest_name = src_orig.replace(OriginPath,'');   
  dest_orig = DestPath + dest_name; print("\t[apply_patch] DstOrig:'%s'" % dest_orig )
  dest_path,dest_name = os.path.split(dest_orig)  
  assert( os.path.exists(dest_path) )

  os.symlink(src_orig, dest_orig)
  print("\t[apply_patch] '%s' patched (New) " %(src_name) )
  del OriginFiles[src_name]

## Checking if all files are patched ...
assert( not len(OriginFiles) ) 


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
