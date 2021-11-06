#include<iostream>
#include<math.h>
#include<random>
#include<time.h>
#include<array>

#define listN std::array<short, N>
const int rootN = 50;
const int N = rootN * rootN;

class lattice2d
{
private:
    std::random_device rand_dev;
    typedef std::mt19937 Gen;
    Gen generator;
    std::uniform_int_distribution<int> n_uniform;
    std::uniform_real_distribution<double> r_uniform;

    double J = 1;
    double k_boltzmann = 1;

    double weights[5];

public:
    listN spins;
    double energy;
    double mag_total;
    double temperature;

public:
    lattice2d(double t);
    void setTemperature(double t);
    void printConfiguration();
    void randomResample();
    void setSpins(listN values);
    int randomFlip();
    void reachEquilibrium(int steps);
    int getIndex(int row, int col);
    short getNearestSum(int index);
    double getEnergy();
    double getM();
    double getm();
};