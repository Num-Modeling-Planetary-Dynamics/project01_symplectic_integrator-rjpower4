#include <array>
#include <cmath>

class Vector3 {
public:
  Vector3(double x, double y, double z) : data_{x, y, z} {}

  double X() const;
  double Y() const;
  double Z() const;

  double Norm() const;

  Vector3 Direction() const;

private:
  std::array<double, 3> data_;
};

double Vector3::X() const { return data_[0]; }
double Vector3::Y() const { return data_[1]; }
double Vector3::Z() const { return data_[2]; }

double Vector3::Norm() const {
  double x = X();
  double y = Y();
  double z = Z();

  return std::sqrt(x * x + y * y + z * z);
}

class State {
private:
  Vector3 position_;
  Vector3 velocity_;
};
