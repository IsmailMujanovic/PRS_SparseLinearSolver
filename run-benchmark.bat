@echo off
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run-benchmark.ps1" %*
