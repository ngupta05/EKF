# EKF
Extended Kalman Filter

This repo implements extended Kalman filter for Udacity's Self Driving Car Nanodegree project.

Brief explanation of what each source file does:
- src/main.cpp: Reads in input file line by line, creates Lidar/Radar update and call EKF's processLidar/Radar functions. Saves the
updated state after each step, to compute RMSE later
- EKF.h/cpp: Implements Extended Kalman Filter logic and also does fusion between Radar and Lidar updates
- Matrix.h/cpp: Basic implementation of matrix datastructure
- State.h: Describes state for interaction between EKF and main.cpp
- Tools.h/cpp: Implements RMSE calculator
- Update.h: Defines datastructures for Lidar and Radar measurements
