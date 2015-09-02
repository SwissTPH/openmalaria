param(
[string]$dir,
[System.Uri]$src
)

$install_path = $dir+"\"
$gsl_lib = $install_path+'gsl-libs.zip'
$gsl_dl = $src.AbsoluteUri

if((Test-Path $gsl_lib) -eq $false) {
    Write-Verbose "Downloading compiled gsl library"
    Write-Verbose $gsl_dl' -> '$gsl_lib
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
    Write-Verbose "Extracting files"
    foreach($item in $zip.items())
    {
      $shell.NameSpace($lib).copyhere($item)
    }
    Write-Verbose "Done."
}
