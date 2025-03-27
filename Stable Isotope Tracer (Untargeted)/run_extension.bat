@echo off

echo Activating environment...
call venv\Scripts\activate

echo Running SIRIUS workflow extension...
python SIRIUS_Extension.py

echo.
echo Deactivating environment...
deactivate

echo ==========================================
echo Workflow run is complete!
echo ==========================================
pause
