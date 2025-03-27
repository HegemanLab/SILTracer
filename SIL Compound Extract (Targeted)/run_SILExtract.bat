@echo off

echo Activating environment...
call venv\Scripts\activate

echo Running the workflow...
python Run.py

echo.
echo Deactivating environment...
deactivate

echo ==========================================
echo Workflow run is complete!
echo ==========================================
pause
