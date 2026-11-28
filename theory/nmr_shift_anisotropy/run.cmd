@echo off
setlocal
cd /d "%~dp0"
if not exist target\classes mkdir target\classes
if defined JAVA_HOME (
  set "NMR_JAVA=%JAVA_HOME%\bin\java.exe"
  set "NMR_JAVAC=%JAVA_HOME%\bin\javac.exe"
) else (
  set "NMR_JAVA=java.exe"
  set "NMR_JAVAC=javac.exe"
)
> target\sources.list (
  for %%f in (src\main\java\org\eljah\nmr\*.java) do @echo src/main/java/org/eljah/nmr/%%~nxf
)
"%NMR_JAVA%" -m jdk.compiler/com.sun.tools.javac.Main -encoding UTF-8 --release 17 -d target\classes @target\sources.list
if errorlevel 1 exit /b 1
"%NMR_JAVA%" -Xmx1500m -cp target\classes org.eljah.nmr.NmrShiftAnisotropy %*
exit /b %ERRORLEVEL%
