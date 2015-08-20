@echo on

set install="C:\projects\openmalaria"

echo %install%

title Installing openMalaria dependencies to %install%

REM install dependencies for openMalaria
REM xsd, xerces-c, boost by download
REM gsl and zlib as submodules

REM Installing boost
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\boost.ps1' -dir %install% -src 'https://github.com/boostorg/boost/archive/boost-1.58.0.zip' -version '1_58_0'"
dir %install%\boost*\*
mkdir %install%\boost_1_58_0
move %install%\boost-boost-1.58.0\libs %install%\boost_1_58_0\

mkdir %install%\lib
REM Installing gsl libs
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\gsl-libs.ps1' -dir %install% -src 'https://github.com/tph-thuering/gsl/releases/download/gsl/gsl-libs.zip'"

mkdir %install%\gsl\gsl\
REM Installing gsl headers
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\gsl-headers.ps1' -dir %install% -src 'https://github.com/tph-thuering/gsl/releases/download/gsl/gsl-headers.zip'"

REM Installing zlib
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\zlib.ps1' -dir %install% -src 'https://github.com/tph-thuering/zlib/releases/download/zlib/zlib.lib.zip'"

REM Installing xsd + xerces-c
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& '.\xsd.ps1' -dir %install% -src 'http://www.codesynthesis.com/download/xsd/4.0/windows/i686/xsd-4.0.msi'"


