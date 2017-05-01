#pragma once
#include "Matrix.h"
#include "State.h"
#include "Update.h"

#include <stdint.h>
#include <vector>

class EKF {
  Matrix m_state;
  bool m_firstUpdate;
  float m_lastTimestamp;
  Matrix m_P;
  Matrix m_F;
  Matrix m_Q;
  Matrix m_HRadar;
  Matrix m_HLidar;
  Matrix m_HLidarT;
  Matrix m_RRadar;
  Matrix m_RLidar;
  Matrix m_I;
  int m_noiseAx;
  int m_noiseAy;

  public:
  EKF();
  void processLidarUpdate(const LidarUpdate& update);
  void processRadarUpdate(const RadarUpdate& update);
  State getState();

  private:
  void updateF(float delta);
  void updateQ(float delta);
  void predict();
  void updateHRadar();
};
