from __future__ import division
import os
import sys
import getopt
try:
  import numpy as np
except ImportError:
  print 'Numpy is not installed'
  sys.exit()
try:
  from PIL import Image
except ImportError:
  print 'PIL not installed'
  sys.exit()

inputTypes = [None, 'tif', 'png']
inputBits = [None, '8', '32']


def convert(location, bits, type, recursNo):
  #nArgs = len(argv)
  # if nArgs < 1:
  #     print 'bin2tiff.py <input file/folder> -t \
  #            <output type (tiff/png)> -r <recursion depth>'
  #     return
  #
  # if argv[0] in ('-h', '--help'):
  #     print "Program to convert the \'.bin\' output from FelixSim 1.7 to image format"
  #     print "First or last argument MUST be path or file"
  #     print "Options:"
  #     print "\t --type (-t) specifies output format"
  #     print "\t\t Supports \'tif\' or \'png\' (default \'tif\')"
  #     print "\t --bits (-b) specifies bit depth for tiff files"
  #     print "\t\t Supports \'8\' or \'32\' (default \'32\')"
  #     print "\t --recursion (-r) specifies how deep in the file structure to search"
  #     print "\t\t (default \'0\')"
  #     return

  if os.path.exists(location):
    print("Found input " + location),
    fname = location.rstrip(os.path.sep)

  ftype = None
  if os.path.isdir(fname):
    ftype = 0
    print "as folder"
  elif os.path.splitext(fname)[1] == '.bin':
    ftype = 1
    print "as file"

  if ftype is None:
    print "Cannot interperet file/folder input"
    return

  strRecurs = recursNo
  recurs = None
  outType = type
  strbitDepth = bits

  valid = True

  if outType not in inputTypes:
    valid = False
    print 'Output format ' + outType + ' not recognised'

  if strbitDepth not in inputBits:
    valid = False
    print 'Bit depth ' + strbitDepth + ' not allowed'

  if ftype == 1 and strRecurs is not None:
    print 'Recursion depth not used for single file'

  if strRecurs is None:
    strRecurs = '0'

  if strRecurs.isdigit():
    recurs = int(strRecurs)
    if recurs < 0:
      valid = False
      print 'Recursion depth must be positive'
  else:
    valid = False
    print "Recursion depth must be an integer"

  if not valid:
    return

  if outType is None:
    print 'No output format selected, using tiff'
    outType = 'tif'

  if outType == 'png' and strbitDepth is not None:
    strbitDepth = '8'
    print 'Bit depth is not used for png files, ignoring'
  elif outType == 'png' and strbitDepth is None:
    strbitDepth = '8'
  elif strbitDepth is None:
    strbitDepth = '32'
    print 'No tiff bit depth set, using 32-bit'

  if ftype == 0:
    doFolder(fname, outType, int(strbitDepth), recurs)
  elif ftype == 1:
    toTiff(fname, int(strbitDepth), outType)
    return


def doFolder(rootdir, ftype, bits, rec):
  count = 0
  startdepth = rootdir.count(os.path.sep)
  for root, subFolders, files in os.walk(rootdir):
    if root.count(os.path.sep) - startdepth > rec:
      continue
    else:
      print root
      for f in files:
        if os.path.splitext(f)[-1].lower() == '.bin':
          count += 1
          toTiff(root + os.path.sep + f, bits, ftype)
  print 'Converted ' + str(count) + ' files'


def toTiff(fname, bits, ftype):
  sz = os.path.getsize(fname)
  sz = sz / 8  # as inputs are 64-bit
  sz = np.sqrt(sz)

  newName = os.path.splitext(fname)[0] + '.' + ftype

  if not sz.is_integer():
    print fname + ' is not a square image, aborting'
    return

  data = np.fromfile(fname, dtype='float64')
  data.resize(sz, sz)

  if bits == 8:
    output = np.uint8(float2int(data, bits))
  # elif bits == 16:
  #     output = np.uint16(float2int(data, bits))
  elif bits == 32:
    output = data.astype('float32')

  Image.fromarray(output).save(newName)


def float2int(data, bits):  # might be really dodgy?
  data -= np.amin(data, axis=None)
  data = data / (np.amax(data) / (2 ** bits - 1))
  return data
