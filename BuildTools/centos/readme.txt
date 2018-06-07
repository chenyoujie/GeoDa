*****************************************************************
*** Build Instructions for GeoDa.  Current as of GeoDa 1.5.23 ***
*****************************************************************

Overview: We assume the build machine hosts a recently-installed
clean OS.  This build guide contains notes on setting up the compile
environment, obtaining the GeoDa source files and dependent libraries,
compiling libraries and GeoDa, and finally packaging the program
for distribution and installation.

***************************************************
*** Building GeoDa for 64-bit CentOS 6 or later ***
***************************************************

NOTE: This is just basic placeholder for now!  Not currently complete.

Build machine assumptions:
- clean CentOS 64-bit installation with all OS updates

1. Install C++ developer tools along with command-line subversion

2. Use SVN to check out GeoDa trunk:
 - From user's home directory: ~/
 - svn co https://geodacenter.repositoryhosting.com/svn/geodacenter_geoda/trunk trunk
 
3. cd to ~/trunk/BuildTools/centos

4. run ./build64.sh to download and build GeoDa and everything it depends upon

5. Package GeoDa for distribution / installation.

