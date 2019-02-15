#include <thread>
#include <functional>
#include <chrono>
#include <cstring>
#include <algorithm>

#include "CTFunc.h"
#include "MRFunc.h"


using namespace std::literals::complex_literals;


double integral_CT(int i, int j, int detectors, int angles)
{
    int steps = 100;

    int ii1 = std::floor(i/detectors); // Angle
    int ii2 = i % detectors; // Detector
    int jj1 = std::floor(j/detectors); // Angle
    int jj2 = j % detectors; // Detector

    double result = 0.0;
    double interval = 2.0/steps;
    double theta = M_PI/angles;

    for (int ii=0;ii<steps;ii++)
    {
        double x = -1.0+ii*interval;
        for (int jj=0;jj<steps;jj++)
        {
            double y = -1.0+jj*interval;

            result += CircleFunc(x, y, 1.0)*
                      DetectorRect(x, y, ii1*theta, ii2, detectors, angles)*
                      DetectorRect(x, y, jj1*theta, jj2, detectors, angles);
        }
    }

    return result;
}


std::complex<double> integral_MR(int i, int j, int detectors, int angles)
{
    int steps = 1000;

    int ii1 = std::floor(i/detectors); // Angle
    int ii2 = i % detectors; // Detector
    int jj1 = std::floor(j/detectors); // Angle
    int jj2 = j % detectors; // Detector

    std::complex<double> result = std::complex<double>(0.0, 0.0);
    double interval = 1.0/steps;
    double theta = M_PI/angles;

    double vx1 = std::cos((double)theta*ii1);
    double vy1 = std::sin((double)theta*ii1);
    double px1 = (double)(ii2 - detectors/2.0+1.0)*vx1;
    double py1 = (double)(ii2 - detectors/2.0+1.0)*vy1;

    double vx2 = std::cos((double)theta*jj1);
    double vy2 = std::sin((double)theta*jj1);
    double px2 = (double)(jj2 - detectors/2.0+1.0)*vx2;
    double py2 = (double)(jj2 - detectors/2.0+1.0)*vy2;

    for (int ii=0;ii<steps;ii++)
    {
        double x = -0.5+ii*interval;
        for (int jj=0;jj<steps;jj++)
        {
            double y = -0.5+jj*interval;

            result += std::exp(-2.0*M_PI*1i*(px1*x+py1*y))*
                      std::exp(-2.0*M_PI*1i*(px2*x+py2*y))*
                      CircleFunc(x, y, 1.0);
        }
    }

    return result;
}


int main()
{
    int detectors = 32;
    int angles = 16;

    std::vector<std::complex<float>> Gram;
    Gram.resize(detectors*angles * detectors*angles);

    for (int i=0;i<detectors*angles;i++)
    {
        for (int j=0;j<detectors*angles;j++)
        {
            Gram[j+i*detectors*angles] = integral_MR(i, j, detectors, angles);
        }
        std::cout<<"i: "<<i<<std::endl;
    }

    std::vector<int> dim =
    {
        detectors*angles,
        detectors*angles
    };
    SaveAsGadgetronRaw("CT_Gram.cplx", Gram.data(), dim);

    return 0;
}
