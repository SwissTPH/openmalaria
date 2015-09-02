param(
[string]$dir,
[System.Uri]$src
)

$install_path = $dir+"\"
$zlib_lib = $install_path+'zlib.lib.zip'
$zlib_dl = $src.AbsoluteUri

if((Test-Path $zlib_lib) -eq $false) {
    Write-Verbose "Downloading compiled zlib library"
    Write-Verbose $zlib_dl' -> '$zlib_lib
    (new-object System.Net.WebClient).DownloadFile($zlib_dl,$zlib_lib)
}

$lib = $install_path+"lib\"
$shell = new-object -com shell.application
$zip = $shell.NameSpace($zlib_lib)
if((Test-Path $lib) -eq $false) {
    New-Item $lib -type directory
}
$zlib_file = $lib+"zlib.lib"
if((Test-Path $zlib_file) -eq $false) {
    Write-Verbose "Extracting zip archive"
    foreach($item in $zip.items())
    {
      $shell.NameSpace($lib).copyhere($item)
    }
    Write-Verbose "Done."
}
