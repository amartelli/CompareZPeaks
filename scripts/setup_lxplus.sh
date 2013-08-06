source /afs/cern.ch/sw/lcg/contrib/gcc/4.3/x86_64-slc5/setup.sh

export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.00/x86_64-slc5-gcc43-opt/root/
export PATH=${ROOTSYS}/bin:${PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

if [ -n "${COMPAREZPEAKS}" ]; then
echo "already set"
else
export THISDIR=`pwd`

export LD_LIBRARY_PATH=${THISDIR}/../CommonUtils/lib:${LD_LIBRARY_PATH}
#export DYLD_LIBRARY_PATH=${THISDIR}/../CommonUtils/lib:${DYLD_LIBRARY_PATH}
export PATH=${THISDIR}/../CommonUtils/bin:${PATH}

export LD_LIBRARY_PATH=${THISDIR}/lib:${LD_LIBRARY_PATH}
#export DYLD_LIBRARY_PATH=${THISDIR}/lib:${DYLD_LIBRARY_PATH}
export PATH=${THISDIR}/bin:${PATH}

export COMMONUTILS=${THISDIR}/../CommonUtils/
export COMMONUTILSINCLUDE=${THISDIR}/../CommonUtils/interface
export COMMONUTILSLIB=${THISDIR}/../CommonUtils/lib

export COMPAREZPEAKS=${THISDIR}/
export COMPAREZPEAKSINCLUDE=${THISDIR}/interface
export COMPAREZPEAKSLIB=${THISDIR}/lib
fi
