#!/usr/bin/env bash
cd "$(dirname "$0")"

echo "Activating virtual environment..."
source venv/bin/activate

echo "Running the SIRIUS workflow extension..."
python SIRIUS_Extension.py

echo
echo "Deactivating environment..."
deactivate

echo "=========================================="
echo "Workflow Extension run is complete!"
echo "=========================================="

read -p "Press Enter to exit."
