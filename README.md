# Introduction
Some tests for the levmar c library. The levmar library is pulled as a git submodule.

# Pull the repo with submodules
```bash
git clone --recurse-submodules git@github.com:spetroce/levmar_test.git
```
*or*
```bash
git clone git@github.com:spetroce/levmar_test.git
cd levmar_test #ie. change to root repo directory
git submodule update --init --recursive
```
If you do not have a github account setup with an SSH key, which is likely the case for EE and ME team members, replace the URL with https://github.com/spetroce/levmar_test.git. So, the shorthand clone instruction is:
```bash
git clone --recurse-submodules https://github.com/spetroce/levmar_test.git
```

# Build Prerequisites
Prerequisite packages to install to perform any building:

Ubuntu:
```bash
sudo apt-get install cmake build-essential git-lfs
```
Fedora:
```bash
sudo dnf install make automake cmake gcc gcc-c++ git git-lfs libtool
```

# Build Instructions

```bash
cd levmar_test
mkdir build
cd build
cmake ..
make
```
