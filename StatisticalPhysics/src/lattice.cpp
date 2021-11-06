#include"lattice.h"

lattice2d::lattice2d(double t) : generator(rand_dev()), n_uniform(std::uniform_int_distribution<>(0, N - 1)), r_uniform(std::uniform_real_distribution<>(0.0f, 1.0f))
{
    setTemperature(t);
    randomResample();
}

void lattice2d::setTemperature(double t)
{
    temperature = t;

    for (int i = 0; i < 5; i++)
    {
        double dE = -J / 2.0f * (8 * i - 16);
        weights[i] = exp(-dE / k_boltzmann / temperature);
    }
}

void lattice2d::printConfiguration()
{
    char c;
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

void lattice2d::randomResample()
{
    for (int i = 0; i < N; i++)
    {
        if (r_uniform(generator) > 0.5f) spins[i] = 1;
        else spins[i] = -1;
    }
    energy = getEnergy();
    mag_total = getM();
}

void lattice2d::setSpins(listN values)
{
    for (int i = 0; i < N; i++)
    {
        spins[i] = values[i];
    }
    energy = getEnergy();
    mag_total = getM();
}

int lattice2d::randomFlip()
{
    int index = n_uniform(generator);
    double r = r_uniform(generator);

    short new_spin = -spins[index];

    int dSum = new_spin * getNearestSum(index) * 4;
    double dE = -J * 0.5f * (double)dSum;

    double boltzmann_weight = weights[dSum / 8 + 2];

    if (dE <= 0 || r <= boltzmann_weight)
    {
        spins[index] = new_spin;
        mag_total += 2.0f * (double)new_spin;
        energy += dE;
        return 1;
    }
    return 0;
}

void lattice2d::reachEquilibrium(int steps)
{
    for (int i = 0; i < N * steps; i++)
    {
        randomFlip();
    }
    //std::cout << "N=" << N << "   " << "kB*T/J=" << temperature * k_boltzmann / J << std::endl;
    //std::cout << N * steps << " steps has been made to reach thermal equilibrium." << std::endl;
    //std::cout << "e=" << energy / N << "   " << "m=" << (double)mag_total / (double)N << std::endl;
    //std::cout << std::endl;
}

int lattice2d::getIndex(int row, int col)
{
    if (row >= rootN) row = row % rootN;
    if (col >= rootN) col = col % rootN;

    if (row < 0) row = rootN - abs(row) % rootN;
    if (col < 0) col = rootN - abs(col) % rootN;

    return rootN * row + col;
}

short lattice2d::getNearestSum(int index)
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

double lattice2d::getEnergy()
{

    double result = 0;
    for (int i = 0; i < N; i++)
    {
        result += spins[i] * getNearestSum(i);
    }
    return -J * 0.5f * result;
}

double lattice2d::getM()
{
    int result = 0;
    for (int i = 0; i < N; i++)
    {
        result += spins[i];
    }
    return (double)result;
}

double lattice2d::getm()
{
    return mag_total / N;
}