#!/bin/tcsh

echo ''
echo 'Setting up python, ROOT and PyROOT'
echo ''

#set arch = slc5_ia32_gcc434
set arch = slc5_amd64_gcc434
#set roofit_release = 5.27.06

setenv CMS_PATH /uscmst1/prod/sw/cms

#setenv PYTHONDIR /uscmst1/prod/sw/cms/slc5_ia32_gcc434/external/python/2.6.4-cms8
setenv PYTHONDIR /uscmst1/prod/sw/cms/slc5_amd64_gcc434/external/python/2.6.4-cms15

#setenv PATH ${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.28.00b/bin:$PATH
#setenv PATH ${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.30.00/bin:$PATH
setenv PATH ${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.30.02_amd64/bin:$PATH

#setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.28.00b
#setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.30.00
#setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.30.02
setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.30.02_amd64

setenv PYTHONPATH ${ROOTSYS}/lib

setenv LD_LIBRARY_PATH ${PYTHONDIR}/lib:${CMS_PATH}/$arch/external/gcc/4.3.4/lib64:${CMS_PATH}/$arch/external/gcc/4.3.4/lib:${ROOTSYS}:${ROOTSYS}/lib:${LD_LIBRARY_PATH}

setenv ROOT_INCLUDE ${ROOTSYS}/include

#setenv PATH ${PWD}:${PATH}
#chmod u+x exost
