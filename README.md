# Satellite Latency Calculator

This project calculates the latency in communication with a satellite given its orbital parameters. It is designed to help you estimate the time delay for signals traveling between a ground station and a satellite in orbit.

## Methodology/Docs
The calculation method is detailed in the PDF located in the `docs/` folder:

> **[Method to Calculate Latency in Communication with a Satellite Given Orbital Data](/docs/Method%20to%20Calculate%20Latency%20in%20Communication%20with%20a%20Satellite%20Given%20Orbital%20Data.pdf)**

This was my final math project during high school and it explains the algorithms and formulas used by this calculator.

## Features
- Calculates communication latency using satellite orbital parameters
- Accepts orbital data from sources such as NORAD
- Allows customization of ground station latitude and longitude

## Requirements
- Python 3.x
- numpy
- datetime (standard library)

## Installation
Install the required Python package:

```pwsh
pip install numpy
```

## Usage
1. Open `satellite_latency_calculator.py` in your editor.
2. Modify the satellite orbital parameters and the ground station latitude/longitude directly in the code as needed.
3. Run the calculator:

```pwsh
python satellite_latency_calculator.py
```

The output will display the calculated latency and other data based on your input parameters.

## Getting Orbital Parameters
You can obtain satellite orbital parameters from [NORAD](https://celestrak.org/NORAD/elements/) or other satellite tracking sources. Enter these values in the code to perform your calculations.
