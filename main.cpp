#include <iostream>
#include <vector>
#include <algorithm>

#include <math.h>

struct point {
  double x, y;
  point operator + (const point& p) const {
    return {p.x + x, p.y + y};
  }
  point operator - (const point& p) const {
    return {-p.x + x, -p.y + y};
  }
  point operator * (double n) const {
    return {n * x, n * y};
  }
  point& operator += (point& p) {
    x += p.x;
    y += p.y;
    return *this;
  }
  point& operator /= (double n) {
    x /= n;
    y /= n;
    return *this;
  }

  double abs() {
    return fabs(x) + fabs(y);
  }
};

double dist(point a, point b) {
  return (a - b).abs();
}

struct function {
  double f(point p) {
    double x = p.x, y = p.y;
    return sin(y)*exp((1-cos(x))*(1-cos(x))) + 
    cos(x)*exp((1-sin(y))*(1-sin(y))) + (x-y)*(x-y);
  }

  double xmin = 0, xmax = 11;
  double ymin = 0, ymax = 11;
};

struct NelderMeadAlgorithm {
  NelderMeadAlgorithm(double alpha, double beta, double gamma, double sigma) :
    alpha(alpha),
    beta(beta),
    gamma(gamma),
    sigma(sigma) {

  };

  bool stop_condition() {
    return dist(simplex[0], simplex[1]) < 0.00001 && dist(simplex[2], simplex[1]) < 0.00001;
  }

  void apply(point init) {
    prepare(init);
    while (!stop_condition())
    {
      set_order();
      new_point_();
    }
    
    std::cout << simplex[0].x << " " << simplex[0].y << " " << values[0];
  }

  void prepare(point init) {
    simplex.push_back(init);
    simplex.push_back(init + point{1, 0});
    simplex.push_back(init + point{0, 1});
    for (auto p: simplex) {
      values.push_back(f.f(p));
    }
  }

  void set_order() {
    if (simplex.size() != 3) {
      throw std::bad_exception();
    }

    int i_min = 0;
    double f_min = values[0];
    int i_max = 0;
    double f_max = values[0];
    for (int i = 0; i < values.size(); ++i) {
      if (values[i] >= f_max) {
        f_max = values[i];
        i_max = i;
      }
      if (values[i] < f_min) {
        f_min = values[i];
        i_min = i;
      }
    }

    int i_mean = 0;
    for (int i = 0; i < values.size(); ++i) {
      if (i != i_max && i != i_min) {
        i_mean = i;
        break;
      }
    }

    order = std::vector<int>{i_min, i_mean, i_max};
  }

  point reflect(point xl, point y, double a) {
    return xl + (y - xl)*a;
  }

  void new_point_() {
    auto& xh = simplex[order[2]];
    auto& xs = simplex[order[1]];
    auto& xl = simplex[order[0]];

    double& fh = values[order[2]];
    double& fs = values[order[1]];
    double& fl = values[order[0]];

    point c = xl*0.5 + xs*0.5;

    auto xr = reflect(c, xh, -alpha);
    double fr = f.f(xr);

    if (fl <= fr && fr < fs) {
      xh = xr;
      fh = fr;
      return;
    }

    if (fr < fl) {
      auto xe = reflect(c, xr, gamma);
      auto fe = f.f(xe);

      if (fe < fr) {
        xh = xe;
        fh = fe;
      } else {
        xh = xr;
        fh = fr;
      }
      return;
    }
    
    if (fr >= fs) {
      if (fr < fh) {
        auto xc = reflect(c, xr, beta);
        auto fc = f.f(xc);

        if (fc <= fr) {
          xh = xc;
          fh = fc;
          return;
        }
      } else {
        auto xc = reflect(c, xh, beta);
        auto fc = f.f(xc);

        if (fc <= fh) {
          xh = xc;
          fh = fc;
          return;
        }
      }
    }

    xh = reflect(xl, xh, sigma);
    xs = reflect(xl, xs, sigma);
    fh = f.f(xh);
    fs = f.f(xs);
  }

  std::vector<point> simplex;
  std::vector<double> values;
  std::vector<int> order;
  

  function f;

  double alpha;
  double beta;
  double gamma;
  double sigma;
};

int main(int argc, char** argv) { //1 0.5 2 0.5
  double x_init = atof(argv[1]);
  double y_init = atof(argv[2]);
  double alpha = atof(argv[3]);
  double beta = atof(argv[4]);
  double gamma = atof(argv[5]);
  double sigma = atof(argv[6]);

  NelderMeadAlgorithm alg(alpha, beta, gamma, sigma);

  alg.apply({x_init, y_init});

  return 0;
}