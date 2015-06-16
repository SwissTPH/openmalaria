@echo on

set install="C:\projects\openmalaria"

echo %install%

title Installing openMalaria dependencies to %install%

REM install dependencies for openMalaria
REM xsd, xerces-c, boost by download
REM gsl and zlib as submodules

REM Installing gsl
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\gsl.ps1' -dir %install% -src 'https://github.com/tph-thuering/gsl/releases/download/gsl/gsl.lib.zip'"

REM Installing zlib
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\zlib.ps1' -dir %install% -src 'https://github.com/tph-thuering/zlib/releases/download/zlib/zlib.lib.zip'"

REM Installing xsd + xerces-c
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\xsd.ps1' -dir %install% -src 'http://www.codesynthesis.com/download/xsd/4.0/windows/i686/xsd-4.0.msi'"

REM Installing boost
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\boost.ps1' -dir %install% -src 'http://netcologne.dl.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.zip'"

