@echo off
setlocal
cd /d "%~dp0"

set "NQR_JAVA=java"
if defined JAVA_HOME if exist "%JAVA_HOME%\bin\java.exe" set "NQR_JAVA=%JAVA_HOME%\bin\java.exe"
if defined SUPERCELL_JAVA set "NQR_JAVA=%SUPERCELL_JAVA%"
if not exist target\classes mkdir target\classes
if errorlevel 1 exit /b %errorlevel%

if defined SUPERCELL_JAVAC goto compiler_override
if defined JAVA_HOME if exist "%JAVA_HOME%\bin\javac.exe" goto compiler_java_home
where javac >nul 2>nul
if not errorlevel 1 goto compiler_path

rem Fallback for a runtime that has jdk.compiler without the javac launcher.
"%NQR_JAVA%" -m jdk.compiler/com.sun.tools.javac.Main -source 17 -target 17 -encoding UTF-8 -d target\classes src\main\java\org\eljah\nmr\SupercellNqr.java src\main\java\org\eljah\nmr\TemperatureSeries.java
goto compiled

:compiler_override
"%SUPERCELL_JAVAC%" --release 17 -encoding UTF-8 -d target\classes src\main\java\org\eljah\nmr\SupercellNqr.java src\main\java\org\eljah\nmr\TemperatureSeries.java
goto compiled

:compiler_java_home
"%JAVA_HOME%\bin\javac.exe" --release 17 -encoding UTF-8 -d target\classes src\main\java\org\eljah\nmr\SupercellNqr.java src\main\java\org\eljah\nmr\TemperatureSeries.java
goto compiled

:compiler_path
javac --release 17 -encoding UTF-8 -d target\classes src\main\java\org\eljah\nmr\SupercellNqr.java src\main\java\org\eljah\nmr\TemperatureSeries.java

:compiled
if errorlevel 1 exit /b %errorlevel%
"%NQR_JAVA%" -cp target\classes org.eljah.nmr.SupercellNqr %*
if errorlevel 1 exit /b %errorlevel%
rem With no explicit arguments, also process the archived temperature CSV.
if not "%~1"=="" exit /b 0
if not exist inputs\temperature_series.csv exit /b 0
"%NQR_JAVA%" -cp target\classes org.eljah.nmr.TemperatureSeries
exit /b %errorlevel%
