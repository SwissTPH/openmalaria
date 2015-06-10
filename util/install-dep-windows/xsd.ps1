param(
[string]$dir
)

$install_path = $dir+"\"
$xsd_msi = $install_path+'xsd-4.0.msi'
$xsd_dl = 'http://www.codesynthesis.com/download/xsd/4.0/windows/i686/xsd-4.0.msi'

if((Test-Path $xsd_msi) -eq $false) {
echo "Downloading xsd package"
echo $xsd_dl
echo $xsd_msi
(new-object System.Net.WebClient).DownloadFile($xsd_dl,$xsd_msi)
}

# Install xsd+xerces-c package
msiexec /passive BASEDIR=$install_path"xsd" /i $xsd_msi /lp xsd.log
