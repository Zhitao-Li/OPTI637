#include <thread>
#include <functional>
#include <chrono>

#include <cstring>

#include "CTFunc.h"
#include "MRFunc.h"


void Process(int pixels, double pixelSize, double integralStep, int detectors, int views, double deltaAngle,
             int view, int detector, double* gridX, double* gridY, std::complex<double>* H,
             std::function<void(double, double, double, int, int, double, double, int, int, std::complex<double>&)> f)
{
    for(int x=0;x<pixels;x++)
    {
        for(int y=0;y<pixels;y++)
        {
            f(pixelSize, integralStep, deltaAngle, view, detector, gridX[x], gridY[y], detectors, views,
              H[y + x*pixels + detector*pixels*pixels + view*detectors*pixels*pixels]);
        }
    }
}



int main(int argc, char** argv)
{
    bool error = false;
    int Modelity = 0;
    int Basis = 0;

    if(argc < 3)
        error = true;

    if(strcmp(argv[1], "-CT") == 0)
        Modelity = 0;
    else if(strcmp(argv[1], "-MR") == 0)
        Modelity = 1;
    else
        error = true;

    if(strcmp(argv[2], "-Delta") == 0)
        Basis = 0;
    else if(strcmp(argv[2], "-Pixel") == 0)
        Basis = 1;
    else if(strcmp(argv[2], "-Fourier") == 0)
        Basis = 2;
    else
        error = true;

    if(error)
    {
        std::cout<<"ForwardModel -CT/-MR -Delta/-Pixel/-Fourier"<<std::endl;
        return -1;
    }

    int detectors = 32;
    int pixels = 128;
    double pixelSize = 2.0/pixels;
    double integralStep = pixelSize/50.0;
    int views = 16;
    double deltaAngle = (double)M_PI/views;

    double* gridX = new double[pixels];
    double* gridY = new double[pixels];

    if(Modelity == 0)
    {
        for (int i=0;i<pixels;i++)
        {
            gridX[i] = -1.0 + pixelSize*i;
            gridY[i] = -1.0 + pixelSize*i;
        }
    }
    else if(Modelity == 1)
    {
        pixelSize = 1.0/pixels;
        for (int i=0;i<pixels;i++)
        {
            gridX[i] = -0.5 + pixelSize*i;
            gridY[i] = -0.5 + pixelSize*i;
        }
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
            if(Modelity == 0)
            {
                if(Basis == 1)
                    t[l] = std::thread(Process,
                                       pixels, pixelSize, integralStep, detectors, views, deltaAngle,
                                       a, l,
                                       gridX, gridY,
                                       H,
                                       CTPixelIntegal);
                else if(Basis == 0)
                    t[l] = std::thread(Process,
                                       pixels, pixelSize, integralStep, detectors, views, deltaAngle,
                                       a, l,
                                       gridX, gridY,
                                       H,
                                       CTDeltaIntegal);
                else if(Basis == 2)
                    t[l] = std::thread(Process,
                                       pixels, pixelSize, integralStep, detectors, views, deltaAngle,
                                       a, l,
                                       gridX, gridY,
                                       H,
                                       CTFourierIntegal);
            }
            else if(Modelity == 1)
            {
                if(Basis == 0)
                    t[l] = std::thread(Process,
                                       pixels, pixelSize, integralStep, detectors, views, deltaAngle,
                                       a, l,
                                       gridX, gridY,
                                       H,
                                       MRDeltaIntegal);
                else if(Basis == 1)
                    t[l] = std::thread(Process,
                                       pixels, pixelSize, integralStep, detectors, views, deltaAngle,
                                       a, l,
                                       gridX, gridY,
                                       H,
                                       MRPixelIntegal);
                else if(Basis == 2)
                    t[l] = std::thread(Process,
                                       pixels, pixelSize, integralStep, detectors, views, deltaAngle,
                                       a, l,
                                       gridX, gridY,
                                       H,
                                       MRFourierIntegal);

            }
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
    if(Modelity == 0)
    {
        if(Basis == 0)
            SaveAsGadgetronRaw("H_CT_Delta.cplx", temp, dim);
        else if(Basis == 1)
            SaveAsGadgetronRaw("H_CT_Pixel.cplx", temp, dim);
        else if(Basis == 2)
            SaveAsGadgetronRaw("H_CT_Fourier.cplx", temp, dim);
    }
    else if(Modelity == 1)
    {
        if(Basis == 0)
            SaveAsGadgetronRaw("H_MR_Delta.cplx", temp, dim);
        else if(Basis == 1)
            SaveAsGadgetronRaw("H_MR_Pixel.cplx", temp, dim);
        else if(Basis == 2)
            SaveAsGadgetronRaw("H_MR_Fourier.cplx", temp, dim);
    }


    delete [] gridX;
    delete [] gridY;
    delete [] H;

    return 0;
}
