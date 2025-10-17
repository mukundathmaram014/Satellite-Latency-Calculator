
# Satellite Latency Calculator

This project calculates the communication latency between a ground station and a satellite using Two-Line Element (TLE) data and user-provided ground coordinates. The C++ program interactively prompts for satellite and ground station information, then computes the round-trip signal latency based on orbital mechanics.

## Methodology/Docs
The calculation method is detailed in the PDF located in the `docs/` folder:

> **[Method to Calculate Latency in Communication with a Satellite Given Orbital Data](/docs/Method%20to%20Calculate%20Latency%20in%20Communication%20with%20a%20Satellite%20Given%20Orbital%20Data.pdf)**

This was my final math project during high school and it explains the algorithms and formulas used by this calculator.


## Features
- Calculates round-trip communication latency using satellite TLE orbital parameters
- Accepts TLE data directly from the user (copy-paste from NORAD or Celestrak)
- Allows entry of ground station latitude and longitude interactively
- Displays detailed satellite and calculation data


## Requirements
- C++17 compatible compiler (e.g., g++, clang++)


## Build Instructions

To compile the C++ program:

```bash
g++ -std=c++17 -o satellite_latency_calculator satellite_latency_calculator.cpp
```

This will produce an executable named `satellite_latency_calculator`.


## Usage
Run the program from your terminal:

```bash
./satellite_latency_calculator
```

You will be prompted to enter:
1. The Satellite TLE data
2. Your current latitude (in degrees)
3. Your current longitude (in degrees)

Example input and output:

```
Enter TLE data for satellite (3 lines, (satellite name on first line and TLE data on next two)):

STARLINK-1012           
1 44718U 19074F   25290.10812837  .00005242  00000+0  37037-3 0  9998
2 44718  53.0517 111.1812 0001607  91.7775 268.3398 15.06419652327103

Enter your current latitude in degrees:
41.878002
Enter your current longitude in degrees:
-93.097702
---------------------------------------
Satellite name: STARLINK-1012           

Inclination: 53.0517
Right Ascension: 111.181
Eccentricity: 0.0001607
Argument of Perigee: 91.7775
Mean Anomaly: 268.34
Mean Motion: 15.0642
Revolutions: 32710
---------------------------------------
Latitude: 42.444
Longitude: -76.5019
---------------------------------------
Would you like to calculate latency for the following data? (Y/n) 
Y

---------------------------------------
Calculated Orbital Data: 

True Anomaly: 113.127°
Satellite-Ground Distance: 5.1019e+06 m
---------------------------------------
Latency: 34.0362 ms
---------------------------------------

```


## Getting Orbital Parameters
You can obtain TLE data for satellites from [NORAD](https://celestrak.org/NORAD/elements/) or other satellite tracking sources. Copy the three TLE lines (name on first line), for your satellite of interest and provide them as input when prompted by the program.

## Accuracy
**Note:** This calculator estimates latency based solely on the physical distance between the ground station and the satellite, using orbital mechanics and the speed of light. It does **not** account for additional sources of delay such as signal processing, atmospheric effects, hardware/software transmission times, or network routing. Actual end-to-end communication latency may be higher in real-world scenarios.

