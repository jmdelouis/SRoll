#!/usr/bin/env python

import os
import sys
import glob
import stools


def usage_and_abort():
  print "usage: %s [--no_rm] <object(s)>" % sys.argv[0]
  print "       --no_rm: only display object names and backend names, don't delete them"
  print "       <object(s)>: full object name(s), accepts wildcards"
  exit(1)
  
if len(sys.argv) == 1:
  usage_and_abort()

params = sys.argv[1:]

remove = True
for i, objname in enumerate( params):
  if objname.strip().lower() == "--no_rm":
    print "--no_rm used, no object will be removed"
    del params[i]
    remove = False
  elif objname.startswith("-"):
    print "error: unkown option " + objname
    print
    usage_and_abort()

for param in params:
  for objname in glob.glob( param):
    print "removing " + objname
    backendname = stools.get_backendname( objname)
    if not os.path.exists( backendname):
      print "error, %s object backendname not found: %s" % (objname, backendname)
    print( "rm -rf " + backendname)
    print( "rm -rf " + objname)
    if remove:
      os.system( "rm -rf " + backendname)
      os.system( "rm -rf " + objname)
