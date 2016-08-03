import sys
import os
import shutil
import subprocess
import glob

# most python users have these; a clear ImportError will be thrown if not
# the original CheckLibraries() function used a 'p_name' variable which was not defined
import numpy
import scipy
try:
    import pyfits
except:
    from astropy.io import fits as pyfits

def getDirs(foldername):
    try:
        return os.walk(foldername).next()[1] # python 2 version
    except:
        return os.walk(foldername).__next__()[1] # python 3 version

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

p = subprocess.Popen('python setup.py build',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
p.wait()
if(p.returncode != 0 and p.returncode != None):
	print("     ----------------------------------------------------------")
	print("     > ERROR: integration.c couldn't be installed.")
	print("     > The error was:\n")
	out, err = p.communicate()
	print((spaced(err,"\t \t")))
	sys.exit()
libfolder = getDirs('build/.')
for name in libfolder:
    if(name[0:3]=='lib'):
        filename = glob.glob('build/'+name+'/*')
        shutil.copy2(filename[0],'.')
shutil.rmtree('build')
print('     >...done!')

p = subprocess.Popen('python setup2.py build',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
p.wait()
if(p.returncode != 0 and p.returncode != None):
	print("     ----------------------------------------------------------")
	print("     > ERROR: Cfunctions.c couldn't be installed.")
	print("     > The error was:\n")
	out, err = p.communicate()
	print((spaced(err,"\t \t")))
	sys.exit()
libfolder = getDirs('build/.')
for name in libfolder:
    if(name[0:3]=='lib'):
        filename = glob.glob('build/'+name+'/*')
        shutil.copy2(filename[0],'.')
shutil.rmtree('build')
print('     >...done!')
