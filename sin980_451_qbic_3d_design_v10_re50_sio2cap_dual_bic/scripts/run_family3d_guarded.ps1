param(
    [string]$Family = "single_air_hole_slab_reference",
    [double]$MinFreeGB = 7.5,
    [int]$Cores = 1,
    [int]$MeshSize = 7,
    [int]$Neigs = 6,
    [int]$PollSeconds = 30,
    [ValidateSet("full", "balanced", "lite")]
    [string]$Profile = "full"
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$Root = Split-Path -Parent $ScriptDir
$LogDir = Join-Path $Root "outputs\family3d\logs"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$safeFamily = $Family -replace "[^A-Za-z0-9_=-]", "_"
$outLog = Join-Path $LogDir "family3d_${safeFamily}_${stamp}.out.log"
$errLog = Join-Path $LogDir "family3d_${safeFamily}_${stamp}.err.log"

$python = "D:\Anaconda\envs\comsol_env\python.exe"
$script = Join-Path $ScriptDir "v10_family3d_pilot.py"
$argsList = @(
    $script,
    "--family", $Family,
    "--cores", "$Cores",
    "--mesh-size", "$MeshSize",
    "--neigs", "$Neigs",
    "--min-free-gb", "$MinFreeGB",
    "--profile", "$Profile"
)

function Get-FreeGB {
    $os = Get-CimInstance Win32_OperatingSystem
    return [double]($os.FreePhysicalMemory / 1MB)
}

Write-Host "[guard] starting family=$Family profile=$Profile cores=$Cores mesh=$MeshSize neigs=$Neigs minFreeGB=$MinFreeGB"
Write-Host "[guard] logs: $outLog"
$proc = Start-Process -FilePath $python -ArgumentList $argsList -WorkingDirectory (Split-Path -Parent $Root) -RedirectStandardOutput $outLog -RedirectStandardError $errLog -PassThru -WindowStyle Hidden

try {
    while (-not $proc.HasExited) {
        Start-Sleep -Seconds $PollSeconds
        $free = Get-FreeGB
        $stampNow = Get-Date -Format "HH:mm:ss"
        Write-Host ("[guard] {0} free={1:N2} GB process={2}" -f $stampNow, $free, $proc.Id)
        if ($free -lt $MinFreeGB) {
            Write-Host ("[guard] free memory {0:N2} GB below threshold {1:N2} GB; stopping process {2}" -f $free, $MinFreeGB, $proc.Id)
            Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
            throw "Memory guard stopped run. Review logs before continuing."
        }
    }
    $proc.Refresh()
    Write-Host "[guard] process exited code=$($proc.ExitCode)"
}
finally {
    if (Test-Path $outLog) {
        Write-Host "[guard] stdout tail:"
        Get-Content $outLog -Tail 40
    }
    if (Test-Path $errLog) {
        $errLines = Get-Content $errLog -Tail 40
        if ($errLines) {
            Write-Host "[guard] stderr tail:"
            $errLines
        }
    }
}
