if (!(-d bin/)) then
mkdir bin
else
echo "bin/ already exists"
endif

if (!(-d lib/)) then
mkdir lib
else
echo "lib/ already exists"
endif

if (!(-d obj/)) then
mkdir obj
else
echo "obj/ already exists"
endif

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
