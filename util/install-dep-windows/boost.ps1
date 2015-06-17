param(
[string]$dir,
[System.Uri]$src
)

$install_path = $dir+"\"
$boost_token = 'boost_1_58_0'
$boost_zip = $install_path+$boost_token+'.zip'
$boost_dl = $src.AbsoluteUri
$boost_install_dir = $install_path + $boost_token

if((Test-Path $boost_zip) -eq $false) {
  echo "Downloading boost package"
  (new-object System.Net.WebClient).DownloadFile($boost_dl,$boost_zip)
}

if((Test-Path $boost_install_dir) -eq $false) {
    # Install boost zipped package
    $shell = new-object -com shell.application
    $zip = $shell.NameSpace($boost_zip)
    foreach($item in $zip.items())
    {
     $shell.NameSpace($install_path).copyhere($item)
    }
}
