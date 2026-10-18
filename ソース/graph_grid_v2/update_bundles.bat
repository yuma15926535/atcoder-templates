@echo off
rem グラフ・グリッド第二世代テンプレートの同梱部分を一括更新する（ダブルクリックで実行）
rem 確認だけしたいときは、コマンドプロンプトで  update_bundles.bat -Check
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0update_bundles.ps1" %*
pause
