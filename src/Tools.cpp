#include "Tools.h"

#include <iostream>

State Tools::calculateRMSE(std::vector<State> estimates,
    std::vector<State> truth) {
  State rmse;

  if (0 == estimates.size() || estimates.size() != truth.size()) {
    std::cout << "Invalid input" << std::endl;
    return rmse;
  }
  unsigned int size = estimates.size();
  for (unsigned int i = 0; i < size; i++) {
    State diff = estimates[i] - truth[i];
    diff.square();
    rmse += diff;
  }

  rmse.x /= size;
  rmse.y /= size;
  rmse.vx /= size;
  rmse.vy /= size;

  rmse.sqrt();
  return rmse;
}
