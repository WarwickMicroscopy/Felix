#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Felix
#
# Richard Beanland, Keith Evans & Rudolf A Roemer
#
# (C) 2013-17, all rights reserved
#
# Version: :VERSION:
# Date:    :DATE:
# Time:    :TIME:
# Status:  :RLSTATUS:
# Build:   :BUILD:
# Author:  :AUTHOR:
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1. step.
# which platform do we use?
# FELIX comes along with some standard settings for various compilers and
# platforms. Often enough it suffices to choose one from the list below.
# In case you do not find an acceptable configuration, you can change the
# associated "makefile.include*" configuration in "makefiles/"
# The supportiung libraries and some dependencies will be provided in the
# associated sub directory of "lib/"

# GNU-compiler-based options

# 64 BIT gcc/gfortran linux system on INTEL architectures
# PLATFORM=INT64NGNU

# 64 BIT gcc/gfortran linux system on AMD architectures
# PLATFORM=OPT64NGNU


# Intel-compiler-based options

# 64 BIT gcc/gfortran linux system on INTEL architectures
PLATFORM=INT64Nifort

# 64 BIT gcc/gfortran linux system on AMD architectures
# PLATFORM=OPT64Nifort


# PGF-compiler-based options
# 64 BIT pgcc/pgf90 linux system on INTEL architectures
# PLATFORM=INT64Npgf

# 64 BIT pgcc/pgf90 linux system on AMD architectures
# PLATFORM=OPT64Npgf

