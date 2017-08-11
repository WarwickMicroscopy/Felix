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

ifeq ($(MYPLATFORM),INT32GNU)
include $(MYSTARTDIR)/makefiles/makefile.includeINT32GNU

else
ifeq ($(MYPLATFORM),INT64NGNU)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64NGNU

else
ifeq ($(MYPLATFORM),INT64YGNU)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64YGNU

else
ifeq ($(MYPLATFORM),OPT32GNU)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT32GNU

else
ifeq ($(MYPLATFORM),OPT64NGNU)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64NGNU

else
ifeq ($(MYPLATFORM),OPT64YGNU)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64YGNU

else
ifeq ($(MYPLATFORM),ITA64Nifort)
include $(MYSTARTDIR)/makefiles/makefile.includeITA64Nifort

else
ifeq ($(MYPLATFORM),ITA64Yifort)
include $(MYSTARTDIR)/makefiles/makefile.includeITA64Yifort

else
ifeq ($(MYPLATFORM),INT32ifort)
include $(MYSTARTDIR)/makefiles/makefile.includeINT32ifort

else
ifeq ($(MYPLATFORM),INT64Nifort)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Nifort

else
ifeq ($(MYPLATFORM),INT64Yifort)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Yifort

else
ifeq ($(MYPLATFORM),OPT32ifort)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT32ifort

else
ifeq ($(MYPLATFORM),OPT64Nifort)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64Nifort

else
ifeq ($(MYPLATFORM),OPT64Yifort)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64Yifort

else
ifeq ($(MYPLATFORM),INT32pgf)
include $(MYSTARTDIR)/makefiles/makefile.includeINT32pgf

else
ifeq ($(MYPLATFORM),INT64Npgf)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Npgf

else
ifeq ($(MYPLATFORM),INT64Ypgf)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Ypgf

else
ifeq ($(MYPLATFORM),OPT32pgf)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT32pgf

else
ifeq ($(MYPLATFORM),OPT64Npgf)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64Npgf

else
ifeq ($(MYPLATFORM),OPT64Ypgf)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64Ypgf

else
ifeq ($(MYPLATFORM),AIXxlf90)
include $(MYSTARTDIR)/makefiles/makefile.includeAIXxlf90

else
ifeq ($(MYPLATFORM),alpha)
include $(MYSTARTDIR)/makefiles/makefile.includealpha

else
ifeq ($(MYPLATFORM),sun)
include $(MYSTARTDIR)/makefiles/makefile.includesun

else
ifeq ($(MYPLATFORM),hpux)
include $(MYSTARTDIR)/makefiles/makefile.includehpux

else
ifeq ($(MYPLATFORM),irix)
include $(MYSTARTDIR)/makefiles/makefile.includeirix

else
ifeq ($(MYPLATFORM),sgi)
include $(MYSTARTDIR)/makefiles/makefile.includesgi

endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
