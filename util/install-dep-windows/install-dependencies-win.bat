@echo on

set install="C:\om_install"

echo %install%

title Installing openMalaria dependencies to %install%

REM install dependencies for openMalaria
REM xsd, xerces-c, boost by download
REM gsl and zlib as submodules

PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\xsd.ps1' -dir %install%"
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\boost.ps1' -dir %install%"

pause
