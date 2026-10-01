@echo off
rem Builds PHANTAST_Live-1.1.jar with the Java compiler and ImageJ library that come with Fiji.
rem Usage:  build.bat "C:\path\to\Fiji"
setlocal

set "FIJI=%~1"
if "%FIJI%"=="" (
    echo Usage: build.bat "C:\path\to\Fiji"
    exit /b 1
)

set "JAVAC="
for /f "delims=" %%J in ('dir /b /s "%FIJI%\java\javac.exe" 2^>nul') do set "JAVAC=%%J"
if not defined JAVAC (
    where javac >nul 2>nul && set "JAVAC=javac"
)
if not defined JAVAC (
    echo Could not find javac.exe in "%FIJI%\java" or on the PATH.
    exit /b 1
)
for %%J in ("%JAVAC%") do set "JAR=%%~dpJjar.exe"
if "%JAVAC%"=="javac" set "JAR=jar"

set "IJJAR="
for %%I in ("%FIJI%\jars\ij-*.jar") do set "IJJAR=%%~fI"
if not defined IJJAR (
    echo Could not find ij-*.jar in "%FIJI%\jars".
    exit /b 1
)

set "HERE=%~dp0"
set "OUT=%HERE%build"
if exist "%OUT%" rmdir /s /q "%OUT%"
mkdir "%OUT%"

echo Compiling with %JAVAC%
"%JAVAC%" --release 8 -Xlint:-options -cp "%IJJAR%" -d "%OUT%" "%HERE%src\PHANTAST_Live.java" || exit /b 1
copy /y "%HERE%plugins.config" "%OUT%\" >nul
copy /y "%HERE%src\PHANTAST_Live.java" "%OUT%\" >nul
if exist "%HERE%..\LICENSE" copy /y "%HERE%..\LICENSE" "%OUT%\LICENSE.txt" >nul
"%JAR%" cf "%HERE%PHANTAST_Live-1.1.jar" -C "%OUT%" . || exit /b 1

echo Built %HERE%PHANTAST_Live-1.1.jar
echo Copy it into "%FIJI%\plugins" and restart Fiji.
