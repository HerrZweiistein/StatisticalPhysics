#include<iostream>
#include<math.h>
#include<random>
#include<array>

#define listN std::array<short, N>

const int rootN = 50;
const int N = rootN * rootN;
const double J = 1;
const double temperature = 2.26;
const double k_boltzmann = 1;

class lattice2d
{
public:
    listN spins;
    double energy;
    int mag_total;

public:
    lattice2d()
    {
        resample();
    }

    void printConfiguration()
    {   
        char c;
        std::cout << std::endl;
        std::cout << "E=" << energy << std::endl;
        std::cout << "M=" << mag_total << std::endl;
        for (int x = 0; x < rootN; x++)
        {
            for (int y = 0; y < rootN; y++)
            {
                int index = getIndex(x, y);
                
                if (spins[index] == 1) c = '+';
                else c = '-';
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void resample()
    {
        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_real_distribution<double> r_uniform(0.0f, 1.0f);

        for (int i = 0; i < N; i++)
        {
            if (r_uniform(generator) > 0.5f) spins[i] = 1;
            else spins[i] = -1;
        }
        energy = getEnergy();
        mag_total = getM();
    }

    void setSpins(listN values)
    {
        for (int i = 0; i < N; i++)
        {
            spins[i] = values[i];
        }
        energy = getEnergy();
        mag_total = getM();
    }

    void randomFlip()
    {
        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_int_distribution<int> n_uniform(0, N - 1);
        std::uniform_real_distribution<double> r_uniform(0.0f, 1.0f);

        int index = n_uniform(generator);
        double r = r_uniform(generator);

        double old_sum = -2.0f * energy / J;
        short new_spin = -spins[index];

        double new_sum = old_sum + 4.0f * new_spin * getNearestSum(index);
        double new_energy = -J * 0.5f * new_sum;

        double delta_energy = new_energy - energy;

        if (delta_energy <= 0 || r <= exp(-delta_energy / (k_boltzmann * temperature)))
        {
            spins[index] = new_spin;
            mag_total += 2 * new_spin;
            energy = new_energy;
        }
    }

    void reachEquilibrium(int steps)
    {
        for (int i = 0; i < steps; i++)
        {
            randomFlip();
        }
        std::cout << "N=" << N << std::endl;
        std::cout << "kB*T/J=" << temperature*k_boltzmann/J << std::endl;
        std::cout << steps << " steps has been made to reach thermal equilibrium." << std::endl;
        std::cout << "E/N=" << energy / N << std::endl << "M/N=m=" <<  (double)mag_total/(double)N << std::endl;
        std::cout << std::endl;
    }

    int getIndex(int row, int col)
    {
        if (row >= rootN) row = row % rootN;
        if (col >= rootN) col = col % rootN;

        if (row < 0) row = rootN - abs(row) % rootN;
        if (col < 0) col = rootN - abs(col) % rootN;

        return rootN * row + col;
    }

    short getNearestSum(int index)
    {
        int row, col, top, right, bottom, left;
        row = index / rootN;
        col = index % rootN;

        top = getIndex(row - 1, col);
        right = getIndex(row, col + 1);
        bottom = getIndex(row + 1, col);
        left = getIndex(row, col - 1);

        return (spins[top] + spins[right] + spins[bottom] + spins[left]);
    }

    double getEnergy()
    {

        double result = 0;
        for (int i = 0; i < N; i++)
        {
            result += spins[i] * getNearestSum(i);
        }
        return -J * 0.5f * result;
    }

    int getM()
    {
        int result = 0;
        for (int i = 0; i < N; i++)
        {
            result += spins[i];
        }
        return result;
    }

    double getm()
    {
        return (double)mag_total / (double)N;
    }
};

int main() {

    lattice2d lat;

    lat.reachEquilibrium(100000);
    lat.printConfiguration();
}