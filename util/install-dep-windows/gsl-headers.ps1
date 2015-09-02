param(
[string]$dir,
[System.Uri]$src
)

$install_path = $dir+"\"
$gsl_lib = $install_path+'gsl-headers.zip'
$gsl_dl = $src.AbsoluteUri

if((Test-Path $gsl_lib) -eq $false) {
    Write-Verbose "Downloading gsl headers"
    Write-Verbose $gsl_dl' -> '$gsl_lib
    (new-object System.Net.WebClient).DownloadFile($gsl_dl,$gsl_lib)
}

$lib = $install_path+"gsl\gsl\"
$shell = new-object -com shell.application
$zip = $shell.NameSpace($gsl_lib)

Write-Verbose "creating new directory: "$lib
if((Test-Path $lib) -eq $false ){
    New-Item $lib -type directory
}
$gsl_file = $lib+"gsl_version.h"
if((Test-Path $gsl_file) -eq $false) {
    Write-Verbose 'Extracting files from '+$gsl_file
    foreach($item in $zip.items())
    {
      if(($item.Name -imatch ".*.h")){
          $shell.NameSpace($lib).copyhere($item)
      }
    }
    Write-Verbose 'Done.'
}

