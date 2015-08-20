param(
[string]$dir,
[System.Uri]$src
)

$install_path = $dir+"\"
$gsl_lib = $install_path+'gsl-libs.zip'
$gsl_dl = $src.AbsoluteUri

if((Test-Path $gsl_lib) -eq $false) {
    echo "Downloading compiled gsl library"
    echo $gsl_dl' -> '$gsl_lib
    (new-object System.Net.WebClient).DownloadFile($gsl_dl,$gsl_lib)
}

$lib = $install_path+"lib\"
$shell = new-object -com shell.application
$zip = $shell.NameSpace($gsl_lib)
if((Test-Path $lib) -eq $false ){
    New-Item $lib -type directory
}
$gsl_file = $lib+"gsl.lib"
if((Test-Path $gsl_file) -eq $false) {
    foreach($item in $zip.items())
    {
      echo "Copying "$item.Name
      $shell.NameSpace($lib).copyhere($item)
    }
}
