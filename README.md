# Pipeline Installation Guide

This is the installation guide for the PESSTO pipeline, used to reduce the data from the NTT telescope.
This new pipeline is compatible with Python 3 and 64bit archs, and finally solves the issues with the Mac M1 chip.
Thanks to the community that has made IRAF compatible with the latest OS and computers.


___
## IRAF

The installation instructions are taken from here: https://iraf-community.github.io/install.html.
For more information about the latest release at the time I wrote this, check here: https://github.com/iraf-community/iraf/releases.

### System Requirements and Dependencies

The distributed binaries require the readline or libedit, curl, and
expat libraries to be installed.

On Debian and its derivatives (Ubuntu, Mint, Devuan, Raspbian etc.):

```code
sudo apt install gcc make flex
sudo apt install libcurl4-openssl-dev libexpat-dev libreadline-dev
```

On Fedora and its derivatives (Redhat, Scientific Linux etc.)

```code
sudo dnf install gcc make perl flex
sudo dnf install libcurl-devel expat-devel readline-devel
```

On MacOS X, you need to have the XCode tools installed. If you
haven't, you can install them with:

```code
xcode-select --install
```

Click "Install" to download and install Xcode Command Line Tools.

### Compile the Sources

The source distribution file is built as a tarball with the package
name and version as base directory. Thus, distribution files can be
unpacked with the command

```code
tar zxf /<path>/iraf-2.17.tar.gz
cd iraf-2.17/
```

In the source directory, execute the install script to create needed
links:

```code
./install 		# execute the install script
```

The script will prompt you for the path to the default image 
directory, the cache directory and the binary files directory.
Usually, you can everywhere use the default settings when asked from 
the install script. You will need to include the binary files 
directory in your PATH before proceeding to the `<make>` step.
In BASH this can be done with the command:

```code
export PATH=/path/to/iraf/bin/:$PATH
```

where `</path/to/iraf/bin/>` is the binary files path specified to 
the install script.

Now you can configure the system for the proper architecture and build:

```code
make <arch>
make sysgen 2>&1 | tee build.log  # this takes some time (~17 min. for me) and prints lots of warnings
```

For `<arch>`, use the proper IRAF architecture name:

`<arch>`   | Operating system | Supported CPU types
-----------|------------------|---------------------------------------
`linux64`  | Linux 64 bit     | x86_64, arm64, mips64, ppc64, riscv64, alpha
`linux`    | Linux 32 bit     | i386, x32, arm, mips
`macos64`  | macOS 64 bit     | arm64
`macintel` | macOS 64 bit     | x86_64
`macosx`   | macOS 32 bit     | i386
`freebsd64`| FreeBSD 64 bit   | x86_64
`freebsd`  | FreeBSD 32 bit   | i386, arm
`hurd`     | GNU HURD 32 bit  | i386

Note that Cygwin and big endian architectures like macosx/ppc are not
supported anymore.


### Test the Build

IRAF comes with a small set of basic tests to ensure that the build
works fine.  To execute the tests, run:

```code
./test/run_tests
```

The details of the tests are described [here](https://github.com/iraf-community/iraf/blob/main/test/README.md).


___
## Anaconda environment and dependencies

We will create an anaconda environment with python 3, the necessary dependencies taken from the stsci/astroconda channel and the SWarp package:

```code
conda config --add channels http://ssb.stsci.edu/astroconda
conda create -n pessto python=3.7 stsci  # be patient, this takes some time to finish as well
conda install -c conda-forge astromatic-swarp
```

### PyRAF

The installation instructions are taken from here: https://iraf-community.github.io/pyraf.html.
Once the anaconda environment has been created, we can proceed to install PyRAF:

```code
conda activate pessto
pip3 install pyraf
```




___
# PESSTO Pipeline

For now, the best option is to install the pipeline by cloning the repository and using the pessto conda environment with python 3:

```code
git clone https://github.com/svalenti/pessto.git
cd pessto
git checkout python3  # temporary branch for development
cd pessto/trunk
python setup.py install
```

## Test

First, you need to download the test data, which you can do manually from the wiki page or using `gdown` as I show below:

```code
mkdir pessto_test
cd pessto_test
pip install gdown
gdown --id 1KSDqJLKURIoVxvFEfLUPoQmz0x-mHAnv
tar zxf PESSTO_Pipeline_Installation_Test_Data.tgz
rm PESSTO_Pipeline_Installation_Test_Data.tgz
```

Now you can run the test in the usual way:

```code
PESSTOFASTSPEC -i EFOSC.2012-04-12T00\:21\:13.429.fits
```

___
# Common Issues

### Matplotlib backend

**Note that this issue has been fixed in v3.0.0.**
Many **MacOS** users have encountered the same error output (e.g., issues [#46](https://github.com/svalenti/pessto/issues/46) [#52](https://github.com/svalenti/pessto/issues/52), [#53](https://github.com/svalenti/pessto/issues/53), [#57](https://github.com/svalenti/pessto/issues/57)):

```code
...
...
libc++abi.dylib: terminating with uncaught exception of type NSException
Abort trap: 6
PANIC in `/Users/.../noao/bin.macosx/x_apextract.e': Write to IPC with no reader
```

If your error looks similar to this one, make sure that you are using the correct matplotlib backend (**TKAgg**). You can manually add this line every time you import matplotlib:

```code
import matplotlib
matplotlib.use("TKAgg")
```

or modify your `~/.matplotlib/matplotlibrc` file, adding:

```code
backend : TKAgg
``` 

If the file doesn't exist, create one.

___
# Reporting Issues

To report any problem, [open an issue](https://github.com/svalenti/pessto/issues) (preferred option) or contact me directly at t.e.muller-bravo@ice.csic.es.

