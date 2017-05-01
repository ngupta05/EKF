#include "EKF.h"

#include <iostream>
#include <stdexcept>

EKF::EKF():
  m_state(4, 1),
  m_firstUpdate(true),
  m_lastTimestamp(0), m_P(4, 4), m_F(4, 4),
  m_Q(4, 4), m_HRadar(3, 4), m_HLidar(2, 4), 
  m_HLidarT(4, 2), m_RRadar(3, 3), m_RLidar(2, 2), m_I(4,4) {

  m_state(0, 0) = m_state(1, 0) = m_state(2, 0) = m_state(3, 0) = 1;

  // Initialize P
  // Unknown position and speed
  m_P(0, 0) = m_P(1, 1) = 1;
  m_P(2, 2) = m_P(3, 3) = 1000;

  // Initialize F
  m_F(0, 0) = m_F(1, 1) = m_F(2, 2) = m_F(3, 3) = 1;
  // Initialize noiseAx/Ay
  m_noiseAx = m_noiseAy = 9;

  // Initialize Q
  // Q is reset on every update 
  
  // m_HRadar is reset on every update
  // Initialize HLidars
  m_HLidar(0, 0) = m_HLidar(1, 1) = 1;
  m_HLidarT(0, 0) = m_HLidarT(1, 1) = 1;

  // Initialize R matrices
  m_RRadar(0, 0) = m_RRadar(2, 2) = 0.09;
  m_RRadar(1, 1) = 0.0009;
  m_RLidar(0, 0) = m_RLidar(1, 1) = 0.0225;

  // Initialize identity
  m_I(0, 0) = m_I(1, 1) = m_I(2, 2) = m_I(3, 3) = 1;
}

void EKF::updateF(float delta) {
  // Reset apt elements
  m_F(0, 2) = m_F(1, 3) = delta;
}

void EKF::updateQ(float delta) {
  // Reset apt elements
  float dt2 = delta * delta;
  float dt3 = dt2 * delta / 2;
  float dt4 = dt2 * dt2 / 4;
  m_Q(0, 0) = dt4 * m_noiseAx;
  m_Q(0, 2) = dt3 * m_noiseAx;
  m_Q(1, 1) = dt4 * m_noiseAy;
  m_Q(1, 3) = dt3 * m_noiseAy;
  m_Q(2, 0) = dt3 * m_noiseAx;
  m_Q(2, 2) = dt2 * m_noiseAx;
  m_Q(3, 1) = dt3 * m_noiseAy;
  m_Q(3, 3) = dt2 * m_noiseAy;
}

void EKF::updateHRadar() {
  float px = m_state(0, 0);
  float py = m_state(1, 0);
  float vx = m_state(2, 0);
  float vy = m_state(3, 0);

	float d = sqrt(px*px + py*py);

	if (d < 1e-2) {
    return;
	  // throw std::invalid_argument("updateHRadar::Division by zero");
	}

	float d2 = d * d;
	float d3 = d2 * d;
	m_HRadar(0, 0) = px / d;
	m_HRadar(0, 1) = py / d;
	m_HRadar(1, 0) = -py / d2;
	m_HRadar(1, 1) = px / d2;
	m_HRadar(2, 0) = py * (vx*py - vy*px) / d3;
	m_HRadar(2, 1) = px * (vy * px - vx * py) / d3;
	m_HRadar(2, 2) = px / d;
	m_HRadar(2, 3) = py / d;
}

void EKF::predict() {
  m_state = m_F * m_state;
  m_P = m_F * m_P * m_F.transpose() + m_Q;
}

void EKF::processLidarUpdate(const LidarUpdate& update) {
  if (m_firstUpdate) {
    m_firstUpdate = false;
    m_state(0, 0) = update.getX();
    m_state(1, 0) = update.getY();
    m_lastTimestamp = update.getTimestamp();
    return;
  }

  float delta = (update.getTimestamp() - m_lastTimestamp) / 1000000.0;
  m_lastTimestamp = update.getTimestamp();
  updateF(delta);
  updateQ(delta);

  // Predict
  predict();

  // Measure Update
  Matrix z(2, 1);
  z(0, 0) = update.getX();
  z(1, 0) = update.getY();
  
  Matrix y = z - m_HLidar * m_state;
  Matrix S = m_HLidar * m_P * m_HLidarT + m_RLidar;
  Matrix K = m_P * m_HLidarT * S.inverse();
  m_state = m_state + K * y;
  m_P = (m_I - K * m_HLidar) * m_P;
}

void EKF::processRadarUpdate(const RadarUpdate& update) {
  if (m_firstUpdate) {
    m_firstUpdate = false;
    float rho = update.getRho();
    float phi = update.getPhi();
    m_state(0, 0) = rho * cos(phi); // x
    m_state(1, 0) = rho * sin(phi); // y
    m_lastTimestamp = update.getTimestamp();
    return;
  }

  float dt = (update.getTimestamp() - m_lastTimestamp) / 1000000.0;
  m_lastTimestamp = update.getTimestamp();
  updateF(dt);
  updateQ(dt);

  // Predict
  predict();

  // apply measurement
  Matrix z(3, 1);
  z(0, 0) = update.getRho();
  z(1, 0) = update.getPhi();
  z(2, 0) = update.getRhoDot();

  // calculate y
  Matrix hx(3, 1);
  float px = m_state(0, 0);
  float py = m_state(1, 0);
  float vx = m_state(2, 0);
  float vy = m_state(3, 0);
  hx(0, 0) = sqrt(px*px + py*py);
  hx(1, 0) = atan2(py, px);
  hx(2, 0) = (px*vx + py*vy) / hx(0, 0);
  Matrix y = z - hx;
  // Adjust phi in y
  float pi = 3.141592;
  while (y(1, 0) < -pi || y(1, 0) > pi) {
    if (y(1, 0) < -pi)
      y(1, 0) += 2 * pi;
    else
      y(1, 0) -= 2 * pi;
  }

  // update HRadar
  updateHRadar(); // calculates Jacobian
  Matrix S = m_HRadar * m_P * m_HRadar.transpose() + m_RRadar;
  Matrix K = m_P * m_HRadar.transpose() * S.inverse();
  m_state = m_state + K * y;
  m_P = (m_I - K * m_HRadar) * m_P;
}

State EKF::getState() {
  State s;
  s.x = m_state(0, 0);
  s.y = m_state(1, 0);
  s.vx = m_state(2, 0);
  s.vy = m_state(3, 0);
  return s;
}
