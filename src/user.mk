# 1. step.
# which platform do we use?
# ILUPACK comes along with some standard settings for various compilers and
# platforms. Often enough it suffices to choose one from the list below.
# In case you do not find an acceptable configuration, you can change the
# associated "makefile.include*" configuration in "makefiles/"
# The supportiung libraries and some dependencies will be provided in the
# associated sub directory of "libs/"

# GNU-compiler-based options

# 32 BIT gcc/gfortran linux system on INTEL architectures
# PLATFORM=INT32GNU

# 64 BIT gcc/gfortran linux system on INTEL architectures
# PLATFORM=INT64NGNU

# 64 BIT gcc/gfortran linux system with 64 bit integer on INTEL architectures
# PLATFORM=INT64YGNU

# 32 BIT gcc/gfortran linux system on AMD architectures
# PLATFORM=OPT32GNU

# 64 BIT gcc/gfortran linux system on AMD architectures
 PLATFORM=OPT64NGNU

# 64 BIT gcc/gfortran linux system with 64 bit integer on AMD architectures
# PLATFORM=OPT64YGNU


# Intel-compiler-based options

# 32 BIT icc/ifort linux system on INTEL architectures
# PLATFORM=INT32ifort

# 64 BIT gcc/gfortran linux system on INTEL architectures
#  PLATFORM=INT64Nifort

# 64 BIT gcc/gfortran linux system with 64 bit integer on INTEL architectures
# PLATFORM=INT64Yifort

# 64 BIT gcc/gfortran linux system on AMD architectures
# PLATFORM=OPT64Nifort

# 64 BIT gcc/gfortran linux system with 64 bit integer on AMD architectures
# PLATFORM=OPT64Yifort


# PGF-compiler-based options
# 32 BIT pgcc/pgf90 linux system on INTEL architectures
# PLATFORM=INT32pgf

# 64 BIT pgcc/pgf90 linux system on INTEL architectures
# PLATFORM=INT64Npgf

# 64 BIT pgcc/pgf90 linux system with 64 bit integer on INTEL architectures
# PLATFORM=INT64Ypgf

# 32 BIT pgcc/pgf90 linux system on AMD architectures
# PLATFORM=OPT32pgf

# 64 BIT pgcc/pgf90 linux system on AMD architectures
# PLATFORM=OPT64Npgf

# 64 BIT pgcc/pgf90 linux system with 64 bit integer on AMD architectures
# PLATFORM=OPT64Ypgf
