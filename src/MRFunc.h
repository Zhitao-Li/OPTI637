#ifndef MRFUNC_H
#define MRFUNC_H


#include <iostream>
#include <complex>

#include "basicFunc.h"


using namespace std::literals::complex_literals;


void MRFourierIntegal(double pixelSize, double integralStep,double deltaAngle,
                      int view, int detector, double gridX, double gridY,
                      int detectors, int views,
                      std::complex<double>& H)
{
    double vx = std::cos((double)deltaAngle*view);
    double vy = std::sin((double)deltaAngle*view);
    double px = (double)(detector - detectors/2.0+1.0)*vx;
    double py = (double)(detector - detectors/2.0+1.0)*vy;

    double step = 1.0/128.0;

    for(int i=0;i<128;i++)
    {
        for(int j=0;j<128;j++)
        {
            double point[2];
            point[0] = -0.5+i*step;
            point[1] = -0.5+j*step;

            double func = CircleFunc(gridX, gridY, 1.0);

            H += func*std::exp(-2.0*M_PI*1i*(point[0]*gridX+point[1]*gridY)/2.0)*
                      std::exp(-2.0*M_PI*1i*(gridX*px+gridY*py));
        }
    }
    H /= 4.0;
}


void MRPixelIntegal(double pixelSize, double integralStep,double deltaAngle,
                    int view, int detector, double gridX, double gridY,
                    int detectors, int views,
                    std::complex<double>& H)
{
    double vx = std::cos((double)deltaAngle*view);
    double vy = std::sin((double)deltaAngle*view);
    double px = (double)(detector - detectors/2.0+1.0)*vx;
    double py = (double)(detector - detectors/2.0+1.0)*vy;

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


void MRDeltaIntegal(double pixelSize, double integralStep,double deltaAngle,
                    int view, int detector, double gridX, double gridY,
                    int detectors, int views,
                    std::complex<double>& H)
{
    double vx = std::cos((double)deltaAngle*view);
    double vy = std::sin((double)deltaAngle*view);
    double px = (double)(detector - detectors/2.0+1.0)*vx;
    double py = (double)(detector - detectors/2.0+1.0)*vy;

    double upperLeftX = gridX, upperLeftY = gridY;
    double upperRightX = gridX + pixelSize, upperRightY = gridY;
    double lowerLeftX = gridX, lowerLeftY = gridY + pixelSize;
    double lowerRightX = gridX + pixelSize, lowerRightY = gridY + pixelSize;

    double func1 = CircleFunc(upperLeftX, upperLeftY, 1.0);

    double func2 = CircleFunc(upperRightX, upperRightY, 1.0);

    double func3 = CircleFunc(lowerLeftX, lowerLeftY, 1.0);

    double func4 = CircleFunc(lowerRightX, lowerRightY, 1.0);

    if(func1 || func2 || func3 || func4)
        H = std::exp(-2.0*M_PI*1i*(gridX*px+gridY*py));
}


#endif // MRFUNC_H
