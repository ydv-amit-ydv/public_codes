@echo off
rem Start Home Designer on Windows: http://127.0.0.1:8765
cd /d "%~dp0"
python server.py --open %*
