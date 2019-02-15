#include <thread>
#include <functional>
#include <chrono>
#include <cstring>
#include <algorithm>

#include "CTFunc.h"
#include "MRFunc.h"


void Process1(int pixels, double pixelSize, double integralStep, int detectors, int views, double deltaAngle,
             int view, int detector, double* gridX, double* gridY, std::complex<double>* H,
             std::vector<std::complex<double>> traj,
             std::function<void(double, double, double, int, int, double, double, int, int,
                                std::complex<double>&, std::vector<std::complex<double>> traj)> f)
{
    for(int x=0;x<pixels;x++)
    {
        for(int y=0;y<pixels;y++)
        {
            f(pixelSize, integralStep, deltaAngle, view, detector, gridX[x], gridY[y], detectors, views,
              H[y + x*pixels + detector*pixels*pixels + view*detectors*pixels*pixels], traj);
        }
    }
}


void MRSpiralIntegal(double pixelSize, double integralStep,double deltaAngle,
                     int view, int detector, double gridX, double gridY,
                     int detectors, int views,
                     std::complex<double>& H,
                     std::vector<std::complex<double>> traj)
{
    double px = (double)(std::real(traj[detector+view*detectors]));
    double py = (double)(std::imag(traj[detector+view*detectors]));

    for(int i=0;i<std::floor(pixelSize/integralStep);i++)
    {
        for(int j=0;j<std::floor(pixelSize/integralStep);j++)
        {
            double point[2];
            point[0] = gridX+(i+1)*integralStep;
            point[1] = gridY+(j+1)*integralStep;

            double func = PixelFunc(point[0], point[1], gridX, gridY, pixelSize) *
                          CircleFunc(point[0], point[1], 1.0);

            H += func*integralStep*integralStep*std::exp(-2.0*M_PI*1i*(gridX*px+gridY*py));
        }
    }
    H /= pixelSize*pixelSize;
}



int main(int argc, char** argv)
{
    int detectors = 128;
    int views = 128;

    int pixels = 128;
    double pixelSize = 2.0/pixels;
    double integralStep = pixelSize/50.0;
    double deltaAngle = (double)M_PI/views;

    double* gridX = new double[pixels];
    double* gridY = new double[pixels];

    pixelSize = 1.0/pixels;
    for (int i=0;i<pixels;i++)
    {
        gridX[i] = -0.5 + pixelSize*i;
        gridY[i] = -0.5 + pixelSize*i;
    }

    std::complex<double>* H = new std::complex<double>[views*detectors*pixels*pixels]();

    // Trajectory
    double delta = 0.5*M_PI/detectors;
    double b = 1.0;

    std::vector<std::complex<double>> traj;
    traj.resize(detectors*views);
    for (int j=0;j<views;j++)
    {
        for (int i=0;i<detectors;i++)
        {
            double r = b*i*delta;
            double angle = j*2*M_PI/views;
            double r1 = std::cos(angle);
            double r2 = std::sin(angle);
            double r3 = -r2;
            double r4 = r1;

            double x = r1 * r*std::cos(i*delta) + r3 * r*sin(i*delta);
            double y = r2 * r*std::cos(i*delta) + r4 * r*sin(i*delta);
            traj[i+j*detectors] = std::complex<double>(x, y);
        }
    }
    for (int i=0;i<traj.size();i++)
    {
        traj[i] =  traj[i]/std::abs(traj[traj.size()-1])*std::complex<double>(41.0);
    }

    std::vector<std::complex<float>> t;
    t.resize(traj.size());
    for (int i=0;i<traj.size();i++)
    {
        t[i] = traj[i];
    }
    std::vector<int> d = {512};
    SaveAsGadgetronRaw("Traj.cplx", t.data(), d);

    auto start = std::chrono::steady_clock::now();
    for(int a=0;a<views;a++)
    {
        std::thread* t = new std::thread[detectors];
        for(int l=0;l<detectors;l++)
        {
            t[l] = std::thread(Process1,
                               pixels, pixelSize, integralStep, detectors, views, deltaAngle,
                               a, l,
                               gridX, gridY,
                               H, traj,
                               MRSpiralIntegal);
        }
        for(int l=0;l<detectors;l++)
        {
            t[l].join();
        }
        delete [] t;
        std::cout<<"a: "<<a<<std::endl;
    }
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;

    std::cout<<std::chrono::duration<double, std::milli>(diff).count()<<"ms"<<std::endl;


    std::complex<float>* temp = new std::complex<float>[pixels*pixels*detectors*views];

    for(int i=0;i<pixels*pixels*detectors*views; i++)
    {
        temp[i] = H[i];
    }

    std::vector<int> dim =
    {
        pixels*pixels,
        detectors*views
    };

    SaveAsGadgetronRaw("H_MR_Spiral_128_128.cplx", temp, dim);


    delete [] gridX;
    delete [] gridY;
    delete [] H;

    return 0;
}
