# this code cannot be run directly
# do 'source /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin/clik_profile.sh' from your sh shell or put it in your profile

function addvar () {
local tmp="${!1}" ;
tmp="${tmp//:${2}:/:}" ; tmp="${tmp/#${2}:/}" ; tmp="${tmp/%:${2}/}" ;
export $1="${2}:${tmp}" ;
} 

if [ -z "${PATH}" ]; then 
PATH=/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin
export PATH
else
addvar PATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin
fi
if [ -z "${PYTHONPATH}" ]; then 
PYTHONPATH=/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages
export PYTHONPATH
else
addvar PYTHONPATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages
fi
if [ -z "${DYLD_LIBRARY_PATH}" ]; then 
DYLD_LIBRARY_PATH=
export DYLD_LIBRARY_PATH
else
addvar DYLD_LIBRARY_PATH 
fi
if [ -z "${DYLD_LIBRARY_PATH}" ]; then 
DYLD_LIBRARY_PATH=/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib
export DYLD_LIBRARY_PATH
else
addvar DYLD_LIBRARY_PATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib
fi
if [ -z "${DYLD_LIBRARY_PATH}" ]; then 
DYLD_LIBRARY_PATH=
export DYLD_LIBRARY_PATH
else
addvar DYLD_LIBRARY_PATH 
fi
if [ -z "${DYLD_LIBRARY_PATH}" ]; then 
DYLD_LIBRARY_PATH=/Users/drudd/Work/cosmosis/cosmosis/ups/cfitsio/v3_35_0/Darwin64bit+13-prof/lib
export DYLD_LIBRARY_PATH
else
addvar DYLD_LIBRARY_PATH /Users/drudd/Work/cosmosis/cosmosis/ups/cfitsio/v3_35_0/Darwin64bit+13-prof/lib
fi
CLIK_DATA=/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/share/clik
export CLIK_DATA

CLIK_PLUGIN=basic,ffp6_foreground,pep_cib
export CLIK_PLUGIN

