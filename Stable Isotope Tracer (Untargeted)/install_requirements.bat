@echo off

echo Creating virtual environment "venv"...
py -3.11 -m venv venv

echo Activating environment...
call venv\Scripts\activate

echo Upgrading pip...
pip install --upgrade pip

echo Installing pinned dependencies from requirements.txt...
pip install -r requirements.txt

echo.
echo Deactivating environment...
deactivate

echo ==========================================
echo Installation complete!
echo Next, run "run_workflow.bat" to use the workflow.
echo ==========================================
pause
