@echo off
for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
echo %LATEST_VERSION%
REM echo "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"
REM echo "::set-output name=env_bat::C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

