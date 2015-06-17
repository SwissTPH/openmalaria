param(
[string]$dir,
[System.Uri]$src
)

$install_path = $dir+"\"
$xsd_msi = $install_path+'xsd-4.0.msi'
$xsd_dl = $src.AbsoluteUri

if((Test-Path $xsd_msi) -eq $false) {
    echo "Downloading xsd package"
    echo $xsd_dl
    echo $xsd_msi
    (new-object System.Net.WebClient).DownloadFile($xsd_dl,$xsd_msi)
}

$xsd_i = $install_path+"xsd"
if((Test-Path $xsd_i) -eq $false) {
  # Install xsd+xerces-c package
  msiexec /passive BASEDIR=$xsd_i /i $xsd_msi /lp xsd.log
}