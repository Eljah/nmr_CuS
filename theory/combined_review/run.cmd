@echo off
setlocal
cd /d "%~dp0"
if not exist target\combined-review-1.0.0.jar (
  call mvn -B package
  if errorlevel 1 exit /b 1
)
java -jar target\combined-review-1.0.0.jar %*
exit /b %errorlevel%
