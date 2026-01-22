# Thermal Sensor Testing

This repository contains the thermal sensor testing framework for UAV-based wildfire detection and navigation.

## Overview

The project evaluates thermal sensor performance for fire detection using Monte Carlo simulation methods. Path planning algorithms are tested on realistic terrain to assess navigation efficiency and detection capabilities.

## Features

- **Monte Carlo Framework**: Comprehensive simulation framework for thermal sensor testing
- **Realistic Terrain**: Obstacle map generated from OpenStreetMap (OSM) data, matching the Blender environment used in RGBD camera testing
- **Path Planning Algorithms**: Implementation and comparison of three navigation approaches:
  - **RRT** (Rapidly-exploring Random Tree)
  - **A*** (A-star)
  - **D*** (D-star)


## Obstacle Map

The obstacle map is generated from the same OpenStreetMap terrain data used in the Unreal Engine environment, ensuring consistency across simulation platforms. This allows for direct comparison between RGBD camera-based obstacle avoidance and thermal sensor-based fire detection performance.

## Path Planning Algorithms

### RRT (Rapidly-exploring Random Tree)
- Probabilistic path planning
- Suitable for high-dimensional spaces
- Explores unknown environments efficiently

### A* (A-star)
- Heuristic-based search algorithm
- Guarantees optimal path if one exists
- Uses cost and distance estimates

### D* (D-star)
- Dynamic pathfinding algorithm
- Efficient replanning when obstacles change
- Suitable for real-time navigation

## Requirements

- MATLAB 
- Mapping libraries for OSM data processing


## Usage

1. **Load Obstacle Map**: Import the OSM-based terrain data
2. **Configure Monte Carlo Parameters**: Set simulation parameters (number of runs, probability of fire spread)
3. **Run Path Planning Tests**: Execute algorithms (RRT, A*, D*) on the obstacle map
4. **Analyse Results**: Compare performance metrics across algorithms, graphs provided

## Monte Carlo Simulation

The Monte Carlo framework runs multiple simulation iterations to:
- Evaluate thermal sensor detection reliability
- Test path planning robustness under varying conditions
- Generate statistical performance metrics
- Assess fire detection rates across different scenarios


## Related Projects

- [RGBD Camera Testing](../Wildfire-RGBD-Camera-Test/) - Companion project for depth-based obstacle avoidance
