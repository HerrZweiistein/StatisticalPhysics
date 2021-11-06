#include"lattice.h"
#include<string>
#include<sstream>
#include <iomanip>
#include <fstream>
#include <ios>

void write(double T, int sweeps)
{   
    lattice2d lat(T);
    lat.reachEquilibrium(1000);
    double m, m2, e;

    std::string result = "%d";

    for (int j = 0; j < sweeps; j++)
    {
        for (int i = 0; i < N; i++)
        {
            lat.randomFlip();
        }

        m = lat.mag_total / N;
        m2 = m * m;
        e = lat.energy / N;


    }
}

void write(std::ofstream& stream, double val)
{
    stream << std::fixed << std::setprecision(8) << val;
}


int main()
{
    double Tmin = 0.1, Tmax = 3;
    int steps = 800;
    int sweeps = 800;
    int thermalisation = 10000;
    double step_size = (Tmax - Tmin) / steps;

    std::ofstream mstream_e("e.txt"); mstream_e.close();
    std::ofstream mstream_m("m.txt"); mstream_m.close();
    std::ofstream mstream_m2("m2.txt"); mstream_m2.close();

    std::ofstream stream_e("e.txt", std::ios_base::app | std::ios_base::out);
    std::ofstream stream_m("m.txt", std::ios_base::app | std::ios_base::out);
    std::ofstream stream_m2("m2.txt", std::ios_base::app | std::ios_base::out);

    lattice2d lat(Tmin);
    lat.reachEquilibrium(thermalisation);
    listN last_spins = lat.spins;

    for (int i = 0; i < steps; i++)
    {   
        
        double Tcurrent = Tmin + i * step_size;
        lat.setTemperature(Tcurrent);
        lat.setSpins(last_spins);

        std::cout << Tcurrent << std::endl;

        double m, m2, e;

        write(stream_e, Tcurrent);
        write(stream_m, Tcurrent);
        write(stream_m2, Tcurrent);

        for (int j = 0; j < sweeps; j++)
        {
            for (int k = 0; k < N; k++)
            {
                lat.randomFlip();
            }

            m = lat.mag_total / N;
            m2 = m * m;
            e = lat.energy / N;

            stream_e << ";";
            write(stream_e, e);

            stream_m << ";";
            write(stream_m, m);

            stream_m2 << ";";
            write(stream_m2, m2);
            
        }
        stream_e << std::endl;
        stream_m << std::endl;
        stream_m2 << std::endl;

        last_spins = lat.spins;
        //delete &lat;
    }

    stream_e.close();
    stream_m.close();
    stream_m2.close();

}