# Pipeline Installation Guide

This is the installation guide for the ePESSTO+ data-reduction pipeline, used to reduce data (photometry and spectroscopy) obtained with the New Technology Telescope (NTT).
This new pipeline is compatible with Python 3 and 64bit archs, and finally solves the issues with the Mac M1 chip.
Thanks to the community that has made IRAF compatible with the latest OS and computers.

**Thanks to Ken for pointing out that Iraf/Pyraf versions compatible with Python 3 existed. Thanks to Carys and Cosimo that helped test the new version of the pipeline and find solutions to different issues that were found.**

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
The iraf command shortcut also needs to be added.
In BASH this can be done with the command:

```code
export PATH=/path/to/iraf/bin/:$PATH
export iraf=/path/to/iraf/
```

where `</path/to/iraf/bin/>` is the binary files path specified to 
the install script and `</path/to/iraf/>` where iraf is installed.

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


## Anaconda environment and dependencies

We will create an anaconda environment with python 3, the necessary dependencies taken from the stsci/astroconda channel and the SWarp package:

```code
conda config --add channels http://ssb.stsci.edu/astroconda
conda create -n pessto python=3.7 stsci  # be patient, this takes some time to finish as well
conda activate pessto
conda install -c conda-forge astromatic-swarp
pip install PyObjC  # necessary for macOS only
```


## PyRAF

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
cd pessto/trunk
conda activate pessto  # unless you are already using the pessto environment
python setup.py install
```


## Fix cursor issue for macOS

The latest Mac computers with the M1 chip have an issue when using the pyraf display: the cursor freezer when hovered over the display. 
To solve this, I found a workaround in https://github.com/iraf-community/pyraf/issues/107. To make life easier for the user, 
there is a script included in the repository (`fix_cursor_macos.py`). Simply run this script using your anaconda environment used to 
install the pipeline:

```code
conda activate pessto
pip install certifi  # to install certificates for macOS
python fix_cursor_macos.py
```

This will download the `Ptkplot.py` file from the repository and replace your local copy of this file (in your anaconda environment), which 
is the one "causing" the issue. The changes replace the red cross that appears on the pyraf display with a more modest one (a small price for a
solution).

**Note:** If you get any error, please check the [Common Issues](#common-issues) section below.

## SExtractor

If you want to reduce SOFI photometry, [SExtractor](https://www.astromatic.net/software/sextractor/) needs to be installed. 
For MacOS, this can be easily done with the following command (thanks Lluís):

```code
brew install sextractor
```

For Linux systems, this might take a few more steps. I tried `sudo apt-get install sextractor`, but this does not work 
(on Ubuntu 22.04 at least). However, in https://sextractor.readthedocs.io/en/latest/Installing.html is explained how to 
install from source. The steps are summarised below:

```code
# dependencies
sudo apt-get update
sudo apt-get install -y libgl1-mesa-glx sextractor scamp libatlas-base-dev libatlas3-base libfftw3-3 libfftw3-dev libtool autoconf

git clone https://github.com/astromatic/sextractor.git
cd sextractor
sh ./autogen.sh
./configure
make -j
sudo make install
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

Note that this only tests the EFOSC2 reduction part, not SOFI!

___
# Common Issues

**Note that some of these issues should have already been fixed in v3.0.0.**

## Matplotlib backend

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



## Could not import aqutil

**MacOS** needs an additional package which can be installed with the following command:

```code
pip install PyObjC
```

For more information about **PyObjC**, check [this link](https://pyobjc.readthedocs.io/en/latest/).



## CERTIFICATE_VERIFY_FAILED

**MacOS** users might get the following error:

```code
urllib.error.URLError: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate
```

To solve this problem, the certificates need to be installed. Check [this link](https://exerror.com/urllib-error-urlerror-urlopen-error-ssl-certificate_verify_failed-certificate-verify-failed-unable-to-get-local-issuer-certificate/) for possible solutions.


## Unsupported pickle protocol: 5

If you get an error like:

```code
Traceback (most recent call last):
  File "xxxxxxx.py", line x, in <module>
    import pyraf
  ...
  ...
  ...
  File ".../anaconda3/envs/pessto/lib/python3.7/site-packages/pyraf/sqliteshelve.py", line 108, in __getitem__
    return pickle.loads(result[0])
ValueError: unsupported pickle protocol: 5
```

To make life easier for the user, there is a script included in the repository (`fix_pickle_macos.py`) to fix this. Simply run these commands:

```code
conda activate pessto
pip3 install pickle5
python fix_pickle_macos.py
```

or you can manually fix this by going to the file `sqliteshelve.py` (under `~/anaconda3/envs/pessto/lib/python3.7/site-packages/pyraf/`) and change `import pickle` for `import pickle5 as pickle` (hopefully, this should be fixed in future pyraf versions).


## Other Issues

The pipeline installation does not work with miniconda, so full anaconda should be installed.

The default shell in Ventura MacOS is `zsh`, but the pipeline does not seem to work with it. Try switching e.g. to `bash`.

___
# Reporting Issues

To report any problem, [open an issue](https://github.com/svalenti/pessto/issues) (preferred option) or contact me directly at t.e.muller-bravo@ice.csic.es or via Slack.

