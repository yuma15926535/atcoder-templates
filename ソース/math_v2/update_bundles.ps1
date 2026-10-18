<#
=============================================================================
  update_bundles.ps1  ―  グラフ・グリッド第二世代テンプレートの「同梱部分」を一括で更新する
-----------------------------------------------------------------------------
  しくみ
    各テンプレートの本体は、次のガードで囲まれている。
        #ifndef GG_GRAPH_BASE
        #define GG_GRAPH_BASE
          ...
        #endif  // GG_GRAPH_BASE
    ガード名 GG_XXX の「持ち主」は、ファイル名が「数字_xxx.cpp」のファイル
    （例: GG_GRAPH_BASE → 00_graph_base.cpp、GG_ROOTED_TREE → 05_rooted_tree.cpp）。
    持ち主のブロックを正として、ほかのファイルに同梱されている同じガードのブロックを置き換える。

  使い方（このファイルと同じフォルダにテンプレートの .cpp を置く）
    update_bundles.bat をダブルクリック          … 同梱部分を更新する
    update_bundles.bat -Check                     … 確認だけ（ファイルは書き換えない）
    PowerShell から:  .\update_bundles.ps1 [-Check] [-Folder <フォルダ>]

  基盤（例: 00_graph_base.cpp の本体）を直したら、このスクリプトを実行すれば、
  基盤を同梱しているすべてのファイル（05_rooted_tree.cpp など）が同じ内容になる。
  ファイルは UTF-8（BOM 付き）で保存し、元の改行コード（LF / CRLF）を保つ。
=============================================================================
#>
param(
    [switch]$Check,
    [string]$Folder = $PSScriptRoot
)
$ErrorActionPreference = 'Stop'

$utf8bom = New-Object System.Text.UTF8Encoding($true)
$blockRe = New-Object System.Text.RegularExpressions.Regex(
    '(?ms)^#ifndef (GG_[A-Z0-9_]+)[ \t]*\r?\n#define \1[ \t]*\r?\n.*?^#endif[ \t]*//[ \t]*\1[ \t]*(?=\r?\n|\z)')

function Normalize([string]$s) { return $s.Replace("`r`n", "`n") }

$files = @(Get-ChildItem -LiteralPath $Folder -Filter '*.cpp' | Sort-Object Name)
if ($files.Count -eq 0) {
    Write-Host "フォルダに .cpp がありません: $Folder"
    exit 0
}

# ---- 1. 持ち主のブロックを集める ----
$owners = @{}
foreach ($f in $files) {
    $nm = [regex]::Match($f.Name, '^\d+_(.+)\.cpp$')
    if (-not $nm.Success) { continue }
    $guard = 'GG_' + $nm.Groups[1].Value.ToUpper()
    $text = [System.IO.File]::ReadAllText($f.FullName, [System.Text.Encoding]::UTF8)
    foreach ($m in $blockRe.Matches($text)) {
        if ($m.Groups[1].Value -eq $guard) {
            if ($owners.ContainsKey($guard)) {
                Write-Host "[エラー] $guard のブロックが $($f.Name) に 2 つ以上あります。"
                exit 2
            }
            $owners[$guard] = @{ File = $f.Name; Block = (Normalize $m.Value) }
        }
    }
}
Write-Host ("持ち主のブロック: " + (($owners.Keys | Sort-Object) -join ', '))

# ---- 2. 各ファイルの同梱ブロックを持ち主と比べて置き換える ----
$changedFiles = 0
$diffCount = 0
$noOwner = @{}
foreach ($f in $files) {
    $text = [System.IO.File]::ReadAllText($f.FullName, [System.Text.Encoding]::UTF8)
    $eol = "`n"
    if ($text.Contains("`r`n")) { $eol = "`r`n" }
    $sb = New-Object System.Text.StringBuilder
    $last = 0
    $changed = @()
    foreach ($m in $blockRe.Matches($text)) {
        [void]$sb.Append($text, $last, $m.Index - $last)
        $g = $m.Groups[1].Value
        $rep = $m.Value
        if (-not $owners.ContainsKey($g)) {
            $noOwner[$g] = $true
        } elseif ($owners[$g].File -ne $f.Name) {
            if ((Normalize $m.Value) -ne $owners[$g].Block) {
                $changed += $g
                $rep = $owners[$g].Block
                if ($eol -eq "`r`n") { $rep = $rep.Replace("`n", "`r`n") }
            }
        }
        [void]$sb.Append($rep)
        $last = $m.Index + $m.Length
    }
    [void]$sb.Append($text, $last, $text.Length - $last)
    if ($changed.Count -gt 0) {
        $diffCount += $changed.Count
        if ($Check) {
            Write-Host ("[差分あり] " + $f.Name + " : " + ($changed -join ', '))
        } else {
            [System.IO.File]::WriteAllText($f.FullName, $sb.ToString(), $utf8bom)
            Write-Host ("[更新]     " + $f.Name + " : " + ($changed -join ', '))
            $changedFiles++
        }
    }
}
foreach ($g in ($noOwner.Keys | Sort-Object)) {
    Write-Host "[注意] $g の持ち主のファイル（数字_名前.cpp）がフォルダにないので、そのままにしました。"
}
if ($Check) {
    if ($diffCount -eq 0) { Write-Host "すべての同梱部分が持ち主と一致しています。" ; exit 0 }
    Write-Host "差分のあるブロック: $diffCount 個（-Check なしで実行すると更新します）"
    exit 1
}
if ($changedFiles -eq 0) { Write-Host "更新の必要はありませんでした（すべて一致）。" }
else { Write-Host "$changedFiles 個のファイルを更新しました。" }
exit 0
