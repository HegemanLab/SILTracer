#!/usr/bin/env bash
cd "$(dirname "$0")"

echo "Creating virtual environment 'venv'..."
python3.11 -m venv venv

echo "Activating environment..."
source venv/bin/activate

echo "Upgrading pip..."
pip install --upgrade pip

echo "Installing pinned dependencies from requirements.txt..."
pip install -r requirements.txt

echo
echo "Deactivating environment..."
deactivate

echo "=========================================="
echo "Installation complete!"
echo "Next, run './run_workflow_mac.command' to use the workflow."
echo "=========================================="
read -p "Press Enter to exit."
