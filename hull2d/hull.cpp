// Test program for ConvexHull2d class
//
// Command line compile options used:
// Visual Studio 2013: 
//    cl /W4 /O2 /EHsc /Fehull.exe main.cpp
// Linux (Debian) g++ 4.9.2:
//    g++ -std=c++11 -Wall -Wextra -pedantic -O2 -o hull main.cpp
//
// Usage: 
//    hull [COUNT] [TRUE]
//
// Options:
//    COUNT = number of random points to generate for each test. Defaults to 10
//    TRUE =  Flag to write data to file and create gnuplot script. You must specify a COUNT to use this flag.
//            This will create 4 files: dbl_test.dat, dbl_test.gp, int_test.dat, int_test.gp.
// Examples:
//    hull  - Generates 10 random integers and doubles with no write to disk.
//    hull 20  - Generates 20 random integers and doubles with no write to disk.
//    hull 20 true  - Generates 20 random integers and doubles and creates data file and gnuplot script.
//

#include "convexhull2d.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <fstream>
#include <chrono>
#include <string>
#include <limits>



using namespace hpm;
using namespace Eigen;
using Timer = std::chrono::steady_clock;
using Duration = std::chrono::duration<double, std::milli>;
using DoubleLimit = std::numeric_limits<double>;
using IntLimit = std::numeric_limits<int>;

template <typename T>
struct Point2d
{
  T x;
  T y;
};

template <typename T>
void writeData(const std::string& data_file, const std::vector<T>& points, 
               ConvexHull2d<T>& hull, std::vector<std::pair<double, double>> rotated)
{
  std::ofstream data(data_file);
  if (!data.is_open()) return;
  data << std::setprecision(15);
  data << "# Point set\n";
  for (auto i : points) {
    data << i.x() << " " << i.y() << '\n';
  }

  data << "\n# Centroid\n";
  data << hull.getCentroid().first << ' ' << hull.getCentroid().second << "\n";

  data << "\n# Convex hull\n";
  if (!hull.isEmpty()) {
    for (auto i : hull.getHull()) {
      data << i.x() << " " << i.y() << '\n';
    }
    data << hull.getHull().front().x() << " " << hull.getHull().front().y() << '\n';
  }

  data << "\n# Rotated vertices\n";
  if (!rotated.empty()) {
    for (auto i : rotated) {
      data << i.first << " " << i.second << '\n';
    }
    data << rotated.front().first << " " << rotated.front().second << '\n';
  }

}

void writeGPScript(const std::string& script, const std::string& data)
{
  std::ofstream gp(script);
  if (!gp.is_open()) return;
  gp << "set title \"Convex Hull Test\"\n";
  gp << "set key outside right box\n";
  gp << "set size ratio -1\n";
  gp << "plot \"" << data << "\" every :::0::0 title \"All Points\" with points, \\\n";
  gp << "\"\" every :::1::1 title \"Centroid\" with points, \\\n";
  gp << "\"\" every :::2::2 title \"Convex Hull\" with lp, \\\n";
  gp << "\"\" every :::3::3 title \"Rotated\" with lp\n";
  gp << "\npause -1\n";
}

int main(int argc, char* argv[])
{
  size_t num_points = 10; // Default number of points
  bool do_write = false; 

  if (argc > 1) {
    int arg1 = std::strtol(argv[1], nullptr, 10);
    if (arg1 >= 0 && arg1 < IntLimit::max() / 2) {
      num_points = arg1;
    }
    else
      std::cout << "\nInvalid number of points! Defaulting to " << num_points << " points per test.\n";
    if (argc > 2) {
      std::string arg2(argv[2]);
      std::transform(arg2.begin(), arg2.end(), arg2.begin(), ::toupper);
      if (arg2 == "TRUE") do_write = true;
    }
  }

  const int min_int = 0;
  //const int min_int = IntLimit::min();
  const int max_int = IntLimit::max();
  //const double min_double = 0.0;
  const double min_double = -std::sqrt(std::sqrt(DoubleLimit::max())) * 10e24;
  const double max_double = std::sqrt(std::sqrt(DoubleLimit::max())) * 10e24;

  std::mt19937 engine;
  engine.seed(std::random_device{}());

  // Start timer and generate random integers to test
  auto timer_start = Timer::now();
  std::vector<Vector2i> int_points(num_points);
  std::uniform_int_distribution<> int_distribution(min_int, max_int);
  for (auto& i : int_points) {
    i.x() = int_distribution(engine);
    i.y() = int_distribution(engine);
  }

  // Generate random doubles to test
  std::vector<Vector2d> double_points(num_points);
  std::uniform_real_distribution<> double_distribution(min_double, max_double);
  for (auto& i : double_points) {
    i.x() = double_distribution(engine);
    i.y() = double_distribution(engine);
  }
  auto random_time = Timer::now() - timer_start;

  // Create hulls
  auto hull_timer_start = Timer::now();
  ConvexHull2d<Vector2i> int_hull(int_points);
  auto int_hull_time = Timer::now() - hull_timer_start;

  ConvexHull2d<Vector2d> double_hull(double_points);
  auto double_hull_time = Timer::now() - hull_timer_start - int_hull_time;
  auto hull_time = Timer::now() - hull_timer_start;
  
  // Test rotate
  const double to_radians = acos(-1) / 180.0; 
  std::pair<double, double> centroid = double_hull.getCentroid();
  std::vector<std::pair<double, double>> double_rotated = double_hull.rotateVertices(-45.0 * to_radians, centroid);
  std::vector<std::pair<double, double>> int_rotated = int_hull.rotateVertices(45.0 * to_radians);

  if (do_write) {
    writeData("int_test.dat", int_points, int_hull, int_rotated);
    writeGPScript("int_test.gp", "int_test.dat");
    writeData("dbl_test.dat", double_points, double_hull, double_rotated);
    writeGPScript("dbl_test.gp", "dbl_test.dat");
  }
  auto total_time = Timer::now() - timer_start;
  auto write_time = total_time - hull_time - random_time;

  std::cout << "\nMin int: " << min_int << '\n';
  std::cout << "Max int: " << max_int << '\n';
  std::cout << "Min double: " << min_double << '\n';
  std::cout << "Max double: " << max_double << '\n';

  std::cout << "\nRandom creation time: " << Duration(random_time).count() << " ms\n";
  std::cout << "Hull creation time: " << Duration(hull_time).count() << " ms\n";
  std::cout << "Write time: " << Duration(write_time).count() << " ms\n";

  std::cout << "\nTotal time: " << Duration(total_time).count() << " ms\n";
  std::cout << "Total points: " << int_points.size() + double_points.size() << '\n';

  std::cout << "\nInteger hull point count: " << int_hull.getHull().size() << '\n';
  std::cout << "Integer hull creation time: " << Duration(int_hull_time).count() << " ms\n";
  std::cout << "Integer hull area: " << int_hull.getArea() << '\n';
  std::cout << "Integer hull is " << (int_hull.isValid() ? "valid\n" : "not valid\n");

  std::cout << "\nDouble hull point count: " << double_hull.getHull().size() << '\n';
  std::cout << "Double hull creation time: " << Duration(double_hull_time).count() << " ms\n";
  std::cout << "Double hull area: " << double_hull.getArea() << '\n';
  std::cout << "Double hull is " << (double_hull.isValid() ? "valid\n" : "not valid\n");

  std::getchar();
}

