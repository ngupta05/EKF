#pragma once
#include "State.h"

#include <vector>

class Tools {
  public:
    static State calculateRMSE(std::vector<State> estimates,
        std::vector<State> truth);
};
