#include <chrono>
#include <cmath>
#include <iostream>
#include <math.h>
#include <matplot/matplot.h>
#include <random>
#include <thread>
#include <vector>

struct GWOResult {
  double alphaFitness;
  std::vector<double> alphaPos;
  std::vector<double> convergence_curve;

  std::vector<std::vector<std::vector<double>>> history;
};

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0.0, 1.0);

std::vector<std::vector<double>> initialize_positions(int PopSize, int D,
                                                      double LB, double UB) {
  std::vector<std::vector<double>> Positions(PopSize, std::vector<double>(D));

  for (int i = 0; i < PopSize; ++i) {
    for (int j = 0; j < D; ++j) {
      Positions[i][j] = dis(gen) * (UB - LB) + LB;
    }
  }
  return Positions;
}

std::vector<std::vector<double>> initialization(int popSize, int dim, double LB,
                                                double UB) {
  int SS_Boundary = 1;

  auto positions = initialize_positions(popSize, dim, LB, UB);

  return positions;
}

std::vector<std::vector<double>> initialization(int popSize, int dim,
                                                const std::vector<double> &LB,
                                                const std::vector<double> &UB) {
  int SS_Boundary = LB.size();

  std::vector<std::vector<double>> Positions(popSize, std::vector<double>(dim));
  for (int i = 0; i < popSize; ++i) {
    for (int j = 0; j < dim; ++j) {
      Positions[i][j] = dis(gen) * (UB[j] - LB[j]) + LB[j];
    }
  }

  return Positions;
}

GWOResult GWO(int popSize, int MaxT, double LB, double UB, int dim,
              double (*f)(std::vector<double> &x)) {

  GWOResult result;

  std::vector<double> AlphaPos(dim, 0.0);
  double AlphaFitness = INFINITY;

  std::vector<double> BetaPos(dim, 0.0);
  double BetaFitness = INFINITY;

  std::vector<double> DeltaPos(dim, 0.0);
  double DeltaFitness = INFINITY;

  std::vector<std::vector<double>> Positions =
      initialization(popSize, dim, UB, LB);

  std::vector<double> Convergence_curve(MaxT, 0.0);

  int itr = 0;

  while (itr < MaxT) {
    for (int i = 0; i < Positions.size(); ++i) {

      for (int j = 0; j < dim; j++) {
        bool BB_UB = Positions[i][j] > UB;
        bool BB_LB = Positions[i][j] < LB;

        if (BB_UB)
          Positions[i][j] = UB;
        else if (BB_LB)
          Positions[i][j] = LB;
      }

      double Fitness = f(Positions[i]);

      if (Fitness < AlphaFitness) {
        AlphaFitness = Fitness;
        AlphaPos = Positions[i];
      }
      if (Fitness > AlphaFitness && Fitness < BetaFitness) {
        BetaFitness = Fitness;
        BetaPos = Positions[i];
      }
      if (Fitness > AlphaFitness && Fitness > BetaFitness &&
          Fitness < DeltaFitness) {
        DeltaFitness = Fitness;
        DeltaPos = Positions[i];
      }
    }

    double a = 2.0 - (double)itr * (2.0 / (double)MaxT); // control parametr

    for (int i = 0; i < popSize; ++i) {
      for (int j = 0; j < dim; ++j) {
        double r1 = dis(gen);
        double r2 = dis(gen);

        double A1 = 2 * a * r1 - a;
        double C1 = 2 * r2;

        double Dalpha = std::abs(C1 * AlphaPos[j] - Positions[i][j]);
        double X1 = AlphaPos[j] - A1 * Dalpha;

        //------------------------------------

        r1 = dis(gen);
        r2 = dis(gen);

        double A2 = 2 * a * r1 - a;
        double C2 = 2 * r2;

        double Dbeta = std::abs(C2 * BetaPos[j] - Positions[i][j]);
        double X2 = BetaPos[j] - A2 * Dbeta;

        //------------------------------------

        r1 = dis(gen);
        r2 = dis(gen);

        double A3 = 2 * a * r1 - a;
        double C3 = 2 * r2;

        double Ddelta = std::abs(C3 * DeltaPos[j] - Positions[i][j]);
        double X3 = DeltaPos[j] - A3 * Ddelta;

        //------------------------------------

        Positions[i][j] = (X1 + X2 + X3) / (3);
      }
    }

    std::vector<std::vector<double>> currentStep = Positions;
    result.history.push_back(currentStep);

    itr += 1;
    Convergence_curve[itr - 1] = AlphaFitness;
  }

  result.alphaFitness = AlphaFitness;
  result.alphaPos = AlphaPos;
  result.convergence_curve = Convergence_curve;

  return result;
}

double f_sphere(std::vector<double> &x) {
  double sum{};

  for (double val : x) {
    sum += val * val;
  }
  return sum;
}

double f_rastrigin(std::vector<double> &x) {
  double A = 10.0;
  int n = x.size();
  double sum = A * n;

  for (double val : x) {
    //  x^2 - A * cos(2 * pi * x)
    sum += (val * val - A * std::cos(2 * M_PI * val));
  }
  return sum;
}

double f_hyperbolic(std::vector<double> &x) {
  double sum = 0;
  for (double val : x) {
    sum -= 1.0 / (std::abs(val) + 1.0);
  }

  return sum;
}

double f_exponential(std::vector<double> &x) {
  double sum = 0;
  for (double val : x) {
    sum -= std::exp(-std::abs(val));
  }
  return sum;
}

double f_harmonic(std::vector<double> &x) {
  double val = x[0];
  return std::pow(val, 3) * std::pow(3 - val, 5) * std::sin(10 * M_PI * val);
}

double f_izom(std::vector<double> &x) {
  double x1 = x[0];
  double x2 = x[1];
  double exponent = -std::pow(x1 - M_PI, 2) - std::pow(x2 - M_PI, 2);
  return -std::cos(x1) * std::cos(x2) * std::exp(exponent);
}

double f_ackley(std::vector<double> &x) {
  double a = 20, b = 0.2, c = 2 * M_PI;
  double sum1 = x[0] * x[0] + x[1] * x[1];
  double sum2 = std::cos(c * x[0]) + std::cos(c * x[1]);
  return -a * std::exp(-b * std::sqrt(0.5 * sum1)) - std::exp(0.5 * sum2) + a +
         std::exp(1);
}

double f_cross_in_tray(std::vector<double> &x) {
  double fact1 = std::sin(x[0]) * std::sin(x[1]);
  double fact2 =
      std::exp(std::abs(100.0 - std::sqrt(x[0] * x[0] + x[1] * x[1]) / M_PI));
  return -0.0001 * std::pow(std::abs(fact1 * fact2) + 1.0, 0.1);
}

double f_eggholder(std::vector<double> &x) {
  double x1 = x[0];
  double x2 = x[1];
  return -(x2 + 47.0) * std::sin(std::sqrt(std::abs(x1 / 2.0 + (x2 + 47.0)))) -
         x1 * std::sin(std::sqrt(std::abs(x1 - (x2 + 47.0))));
}

double f_holder_table(std::vector<double> &x) {
  double exp_term = std::abs(1.0 - std::sqrt(x[0] * x[0] + x[1] * x[1]) / M_PI);
  return -std::abs(std::sin(x[0]) * std::cos(x[1]) * std::exp(exp_term));
}

double f_schaffer1(std::vector<double> &x) {
  double x2_y2 = x[0] * x[0] + x[1] * x[1];
  double num = std::pow(std::sin(x[0] * x[0] - x[1] * x[1]), 2) - 0.5;
  double den = std::pow(1.0 + 0.001 * x2_y2, 2);
  return 0.5 + num / den;
}

double f_schaffer2(std::vector<double> &x) {
  double x1 = x[0];
  double x2 = x[1];

  double x2_plus_y2 = x1 * x1 + x2 * x2;

  double inner_sin = std::sin(std::abs(x1 * x1 - x2 * x2));

  double num = std::pow(std::cos(inner_sin), 2) - 0.5;

  double den = std::pow(1.0 + 0.001 * x2_plus_y2, 2);

  return 0.5 + num / den;
}

double f1(std::vector<double> &x) { return 0.0; }

int main() {

  bool rastringFunction = false;
  bool hyperbolicFunction = false;
  bool exponentialFunction = false;
  bool izomaFunc = false;
  bool ackleyFunction = false;
  bool crossInTray = false;
  bool eggholderFunction = false;
  bool holderTableFunction = false;
  bool schaffer1Function = true;
  bool schaffer2Function = false;

  GWOResult result;
  double lb, ub;
  int popSize = 100, MaxT = 10, dim;
  auto f_ptr = f1;

  // initialization(50, 10, 0.0, 1.0);

  if (rastringFunction) {
    lb = -5.12;
    ub = 5.12;
    dim = 2;
    f_ptr = f_rastrigin;
  } else if (hyperbolicFunction) {
    lb = -10.0;
    ub = 10.0;
    dim = 2;
    f_ptr = f_hyperbolic;
  } else if (exponentialFunction) {
    lb = -1.0;
    ub = 1.0;
    dim = 2;
    f_ptr = f_exponential;
  } else if (izomaFunc) {
    lb = -100.0, ub = 100.0;
    dim = 2;
    f_ptr = f_izom;
  } else if (ackleyFunction) {
    lb = -5.0;
    ub = 5.0;
    dim = 2;
    f_ptr = f_ackley;
  } else if (crossInTray) {
    lb = -10.0;
    ub = 10.0;
    dim = 2;
    f_ptr = f_cross_in_tray;
  } else if (eggholderFunction) {
    lb = -512.0;
    ub = 512.0;
    dim = 2;
    f_ptr = f_eggholder;
  } else if (holderTableFunction) {
    lb = -10.0;
    ub = 10.0;
    dim = 2;
    f_ptr = f_holder_table;
  } else if (schaffer1Function) {
    lb = -100;
    ub = 100;
    dim = 2;
    f_ptr = f_schaffer1;
  } else if (schaffer2Function) {
    lb = -100.0;
    ub = 100.0;
    dim = 2;
    f_ptr = f_schaffer2;
  }

  auto start = std::chrono::high_resolution_clock::now();

  result = GWO(popSize, MaxT, lb, ub, dim, f_ptr);

  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> duration = end - start;

  auto xGrid = matplot::linspace(lb, ub, 50);
  auto yGrid = matplot::linspace(lb, ub, 50);
  auto [X, Y] = matplot::meshgrid(xGrid, yGrid);

  std::vector<std::vector<double>> Z(X.size(),
                                     std::vector<double>(X[0].size()));
  for (size_t i = 0; i < X.size(); ++i) {
    for (size_t j = 0; j < X[i].size(); ++j) {
      std::vector<double> point = {X[i][j], Y[i][j]};
      Z[i][j] = f_ptr(point);
    }
  }

  std::cout << "--------------------------------" << std::endl;
  std::cout << "AlphaFitness: " << result.alphaFitness << std::endl;
  std::cout << "Time: " << duration.count() << " ms" << std::endl;
  std::cout << "--------------------------------" << std::endl;

  auto fig = matplot::gcf();
  fig->size(1920, 1080);
  fig->quiet_mode(false);

  /*for (size_t t = 0; t < result.history.size(); ++t) {

    std::vector<double> x_coords, y_coords;
    for (const auto &wolf : result.history[t]) {
      x_coords.push_back(wolf[0]);
      y_coords.push_back(dim > 1 ? wolf[1] : 0.0);
    }

    auto s = matplot::scatter(x_coords, y_coords);
    s->marker_size(12);
    s->line_width(2.0);

    // matplot::scatter(x_coords, y_coords);
    matplot::xlim({lb, ub});
    matplot::ylim({dim > 1 ? lb : -2.0, dim > 1 ? ub : 2.0});

    s->x_data(x_coords);
    s->y_data(y_coords);

    fig->draw();
    // std::this_thread::sleep_for(std::chrono::milliseconds(500));
  }*/

  for (size_t t = 0; t < result.history.size(); t++) {
    matplot::cla();

    auto c = matplot::contour(X, Y, Z);
    c->line_width(3.0);
    matplot::hold(matplot::on);

    matplot::colormap(matplot::palette::jet());
    matplot::colorbar();

    std::vector<double> x_coords, y_coords;
    for (const auto &wolf : result.history[t]) {
      x_coords.push_back(wolf[0]);
      y_coords.push_back(dim > 1 ? wolf[1] : 0.0);
    }

    auto s = matplot::scatter(x_coords, y_coords);
    s->marker_size(20);
    s->line_width(5.0);

    s->marker_style(matplot::line_spec::marker_style::circle);
    s->marker_face_color("white");

    // matplot::scatter(x_coords, y_coords);
    matplot::xlim({lb, ub});
    matplot::ylim({dim > 1 ? lb : -2.0, dim > 1 ? ub : 2.0});

    s->x_data(x_coords);
    s->y_data(y_coords);

    fig->draw();
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
  }

  matplot::show();
  return 0;
}
