#include "EKF.h"
#include "State.h"
#include "Tools.h"
#include "Update.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

int main(int argc, char** argv) {
  std::string usage = argv[0];
  usage += " input.txt";
  usage += " output.txt";
  if (3 != argc) {
    std::cout << usage << std::endl;
    return 0;
  }

  const char* input = argv[1];
  const char* output = argv[2];

  std::ifstream inputFile(input, std::ifstream::in);
  std::ofstream outputFile(output, std::ofstream::out);

  std::string line;
  std::string type;
  uint64_t timestamp;
  std::vector<State> groundTruth;
  std::vector<State> estimates;
  EKF ekf;
  while (getline(inputFile, line)) {
    std::istringstream iss(line);

    iss >> type;

    float x, y;
    if (type.compare("L") == 0) {
      iss >> x;
      iss >> y;
      iss >> timestamp;
      LidarUpdate update(x, y, timestamp);
      ekf.processLidarUpdate(update);
    } else {
      float rho, phi, rhodot;
      iss >> rho;
      iss >> phi;
      iss >> rhodot;
      iss >> timestamp;
      RadarUpdate update(rho, phi, rhodot, timestamp);
      ekf.processRadarUpdate(update);
      x = rho * cos(phi);
      y = rho * sin(phi);
    }

    State estimate = ekf.getState();

    State truth;
    iss >> truth.x;
    iss >> truth.y;
    iss >> truth.vx;
    iss >> truth.vy;
    groundTruth.push_back(truth);
    estimates.push_back(estimate);

    outputFile << estimate.x << "\t" << estimate.y << "\t";
    outputFile << estimate.vx << "\t" << estimate.vy << "\t";
    outputFile << x << "\t" << y << "\t";
    outputFile << truth.x << "\t" << truth.y << "\t";
    outputFile << truth.vx << "\t" << truth.vy << "\n";
  }

  State rmse = Tools::calculateRMSE(estimates, groundTruth); 

  std::cout << "Accuracy - RMSE:" << std::endl
    << rmse.x << "\t" << rmse.y << "\t"
    << rmse.vx << "\t" << rmse.vy << "\t" << std::endl;

  return 0;
}
