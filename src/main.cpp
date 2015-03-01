#include <iostream>
#include <thread>
#include <functional>
#include <chrono>

#include "CTFunc.h"


void Process(int pixels, double pixelSize, double integralStep, int detectors, double deltaAngle,
             int view, int detector, double* gridX, double* gridY, std::complex<double>* H,
             std::function<void(double, double, double, int, int, double, double, std::complex<double>&)> f)
{
    for(int x=0;x<pixels;x++)
    {
        for(int y=0;y<pixels;y++)
        {
            f(pixelSize, integralStep, deltaAngle, view, detector, gridX[x], gridY[y],
              H[y + x*pixels + detector*pixels*pixels + view*detectors*pixels*pixels]);
        }
    }
}



int main()
{
    int detectors = 32;
    int pixels = 128;
    double pixelSize = 2.0/pixels;
    double integralStep = pixelSize/50.0;
    int views = 16;
    double deltaAngle = (double)M_PI/views;

    double* gridX = new double[pixels];
    double* gridY = new double[pixels];

    for (int i=0;i<pixels;i++)
    {
        gridX[i] = -1.0 + pixelSize*i;
        gridY[i] = -1.0 + pixelSize*i;
    }

    std::complex<double>* H = new std::complex<double>[views*detectors*pixels*pixels];

    for(int i=0;i<views*detectors*pixels*pixels;i++)
    {
        H[i] = 0;
    }


    auto start = std::chrono::steady_clock::now();
    for(int a=0;a<views;a++)
    {
        std::thread t[detectors];
        for(int l=0;l<detectors;l++)
        {
//            t[l] = std::thread(Process,
//                               pixels, pixelSize, integralStep, detectors, deltaAngle,
//                               a, l,
//                               gridX, gridY,
//                               H,
//                               PixelIntegal);
//            t[l] = std::thread(Process,
//                               pixels, pixelSize, integralStep, detectors, deltaAngle,
//                               a, l,
//                               gridX, gridY,
//                               H,
//                               DeltaIntegal);
            t[l] = std::thread(Process,
                               pixels, pixelSize, integralStep, detectors, deltaAngle,
                               a, l,
                               gridX, gridY,
                               H,
                               FourierIntegal);
        }
        for(int l=0;l<detectors;l++)
        {
            t[l].join();
        }
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
    SaveAsGadgetronRaw("H_fourier.cplx", temp, dim);


    delete [] gridX;
    delete [] gridY;
    delete [] H;

    return 0;
}
