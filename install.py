import sys
import os
import shutil
import subprocess
import glob

def getDirs(foldername):
    return os.walk(foldername).next()[1]

def spaced(input,space):
    fixed = False
    i = 0
    input = space+input
    while(not fixed):
        if(input[i:i+1] == '\n'):
           input = input[0:i+1]+space+input[i+1:]
           i = i + len(space)
        i = i + 1
        if(i == len(input)-1):
          fixed = True
    return input

def CheckLibraries():
    print "     ----------------------------------------------------------"
    try:
      import numpy
    except ImportError:
      print "     ----------------------------------------------------------"
      print '     ERROR: '+p_name+' will not be installed in your system because' 
      print '             numpy is not installed in your system.'
      print '             To install it, go to: http://www.numpy.org/\n\n'
      sys.exit(1)
    print "     > Numpy is ok!"
    try:
      import scipy
    except ImportError:
      print "     ----------------------------------------------------------"
      print '     ERROR: '+p_name+' will not be installed in your system because'
      print '             scipy is not installed in your system.'
      print '             To install it, go to: http://www.scipy.org/\n\n'
      sys.exit(1)
    print "     > Scipy is ok!"
    try:
      import astropy 
    except ImportError:
      print "     ----------------------------------------------------------"
      print '     ERROR: '+p_name+' will not be installed in your system because'
      print '            pyfits is not installed in your system.'
      print '            To install it: pip install astropy \n\n'
      sys.exit(1)
    print "     > astropy is ok!"

CheckLibraries()

p = subprocess.Popen('python setup.py build',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
p.wait()
if(p.returncode != 0 and p.returncode != None):
	print "     ----------------------------------------------------------"
	print "     > ERROR: integration.c couldn't be installed."
	print "     > The error was:\n"
	out, err = p.communicate()
	print spaced(err,"\t \t")
	sys.exit()
libfolder = getDirs('build/.')
for name in libfolder:
    if(name[0:3]=='lib'):
        filename = glob.glob('build/'+name+'/*')
        shutil.copy2(filename[0],'.')
shutil.rmtree('build')
print '     >...done!'

p = subprocess.Popen('python setup2.py build',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
p.wait()
if(p.returncode != 0 and p.returncode != None):
	print "     ----------------------------------------------------------"
	print "     > ERROR: Cfunctions.c couldn't be installed."
	print "     > The error was:\n"
	out, err = p.communicate()
	print spaced(err,"\t \t")
	sys.exit()
libfolder = getDirs('build/.')
for name in libfolder:
    if(name[0:3]=='lib'):
        filename = glob.glob('build/'+name+'/*')
        shutil.copy2(filename[0],'.')
shutil.rmtree('build')
print '     >...done!'
