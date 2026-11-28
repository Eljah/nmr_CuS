@echo off
setlocal
cd /d "%~dp0"
if not defined JOINT_JAVA set "JOINT_JAVA=java"
if not exist target\classes mkdir target\classes
where javac >nul 2>nul
if errorlevel 1 (
  "%JOINT_JAVA%" -m jdk.compiler/com.sun.tools.javac.Main -source 17 -target 17 -encoding UTF-8 -d target\classes src\main\java\org\eljah\nmr\joint\Numerics.java src\main\java\org\eljah\nmr\joint\JointTemperature.java
) else (
  javac --release 17 -encoding UTF-8 -d target\classes src\main\java\org\eljah\nmr\joint\Numerics.java src\main\java\org\eljah\nmr\joint\JointTemperature.java
)
if errorlevel 1 exit /b 1
"%JOINT_JAVA%" -cp target\classes org.eljah.nmr.joint.JointTemperature %*
exit /b %errorlevel%
