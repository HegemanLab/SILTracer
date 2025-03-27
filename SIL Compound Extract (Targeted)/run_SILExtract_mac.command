#!/usr/bin/env bash
cd "$(dirname "$0")"

echo "Activating virtual environment..."
source venv/bin/activate

echo "Running workflow..."
python Run.py

echo
echo "Deactivating environment..."
deactivate

echo "=========================================="
echo "Workflow run is complete!"
echo "=========================================="

read -p "Press Enter to exit."
