# this code cannot be run directly
# do 'source /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin/clik_profile.csh' from your csh shell or put it in your profile

 

if !($?PATH) then
setenv PATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin
else
set newvar=$PATH
set newvar=`echo ${newvar} | sed s@:/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin:@:@g`
set newvar=`echo ${newvar} | sed s@:/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin\$@@` 
set newvar=`echo ${newvar} | sed s@^/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin:@@`  
set newvar=/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin:${newvar}                     
setenv PATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin:${newvar} 
endif
if !($?PYTHONPATH) then
setenv PYTHONPATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages
else
set newvar=$PYTHONPATH
set newvar=`echo ${newvar} | sed s@:/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages:@:@g`
set newvar=`echo ${newvar} | sed s@:/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages\$@@` 
set newvar=`echo ${newvar} | sed s@^/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages:@@`  
set newvar=/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages:${newvar}                     
setenv PYTHONPATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages:${newvar} 
endif
if !($?DYLD_LIBRARY_PATH) then
setenv DYLD_LIBRARY_PATH 
else
set newvar=$DYLD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@::@:@g`
set newvar=`echo ${newvar} | sed s@:\$@@` 
set newvar=`echo ${newvar} | sed s@^:@@`  
set newvar=:${newvar}                     
setenv DYLD_LIBRARY_PATH :${newvar} 
endif
if !($?DYLD_LIBRARY_PATH) then
setenv DYLD_LIBRARY_PATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib
else
set newvar=$DYLD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@:/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib:@:@g`
set newvar=`echo ${newvar} | sed s@:/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib\$@@` 
set newvar=`echo ${newvar} | sed s@^/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib:@@`  
set newvar=/Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib:${newvar}                     
setenv DYLD_LIBRARY_PATH /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib:${newvar} 
endif
if !($?DYLD_LIBRARY_PATH) then
setenv DYLD_LIBRARY_PATH 
else
set newvar=$DYLD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@::@:@g`
set newvar=`echo ${newvar} | sed s@:\$@@` 
set newvar=`echo ${newvar} | sed s@^:@@`  
set newvar=:${newvar}                     
setenv DYLD_LIBRARY_PATH :${newvar} 
endif
if !($?DYLD_LIBRARY_PATH) then
setenv DYLD_LIBRARY_PATH /Users/drudd/Work/cosmosis/cosmosis/ups/cfitsio/v3_35_0/Darwin64bit+13-prof/lib
else
set newvar=$DYLD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@:/Users/drudd/Work/cosmosis/cosmosis/ups/cfitsio/v3_35_0/Darwin64bit+13-prof/lib:@:@g`
set newvar=`echo ${newvar} | sed s@:/Users/drudd/Work/cosmosis/cosmosis/ups/cfitsio/v3_35_0/Darwin64bit+13-prof/lib\$@@` 
set newvar=`echo ${newvar} | sed s@^/Users/drudd/Work/cosmosis/cosmosis/ups/cfitsio/v3_35_0/Darwin64bit+13-prof/lib:@@`  
set newvar=:${newvar}                     
setenv DYLD_LIBRARY_PATH :${newvar} 
endif
setenv CLIK_DATA /Users/drudd/Work/cosmosis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/share/clik

setenv CLIK_PLUGIN basic,ffp6_foreground,pep_cib

