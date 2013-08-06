source /afs/cern.ch/sw/lcg/contrib/gcc/4.3/x86_64-slc5/setup.csh

setenv ROOTSYS /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.00/x86_64-slc5-gcc43-opt/root/
setenv PATH ${ROOTSYS}/bin:${PATH}
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}


if (${?COMPAREZPEAKS}) then
echo "already set"
else
setenv THISDIR `pwd`

setenv LD_LIBRARY_PATH ${THISDIR}/../CommonUtils/lib:${LD_LIBRARY_PATH}
#setenv DYLD_LIBRARY_PATH ${THISDIR}/../CommonUtils/lib:${DYLD_LIBRARY_PATH}
setenv PATH ${THISDIR}/../CommonUtils/bin:${PATH}

setenv LD_LIBRARY_PATH ${THISDIR}/lib:${LD_LIBRARY_PATH}
#setenv DYLD_LIBRARY_PATH ${THISDIR}/lib:${DYLD_LIBRARY_PATH}
setenv PATH ${THISDIR}/bin:${PATH}

setenv COMMONUTILS ${THISDIR}/../CommonUtils
setenv COMMONUTILSINCLUDE ${THISDIR}/../CommonUtils/interface
setenv COMMONUTILSLIB ${THISDIR}/../CommonUtils/lib

setenv COMPAREZPEAKS ${THISDIR}/
setenv COMPAREZPEAKSINCLUDE ${THISDIR}/interface
setenv COMPAREZPEAKSLIB ${THISDIR}/lib
endif
