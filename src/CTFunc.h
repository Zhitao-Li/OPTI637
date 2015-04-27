#ifndef CTFUNC_H
#define CTFUNC_H


#include <iostream>
#include <complex>

#include "basicFunc.h"


using namespace std::literals::complex_literals;


auto DetectorRect = [&](double x, double y, double theta, int l, int detectors, int views)
{
    double rvec[2], nvec[2];
    rvec[0] = x;
    rvec[1] = y;
    nvec[0] = -std::sin(theta);
    nvec[1] = std::cos(theta);

    double dot = rvec[0]*nvec[0] + rvec[1]*nvec[1];
    double halfDetectorSize = (double)detectors/2.0;

    return rect((dot-((l-halfDetectorSize)/halfDetectorSize + 1.0/(double)detectors))/(1.0/halfDetectorSize));
};


void CTFourierIntegal(double pixelSize, double integralStep,double deltaAngle,
                      int view, int detector, double gridX, double gridY,
                      int detectors, int views,
                      std::complex<double>& H)
{
    double step = 2.0/128.0;

    for(int i=0;i<128;i++)
    {
        for(int j=0;j<128;j++)
        {
            double point[2];
            point[0] = -1.0+i*step;
            point[1] = -1.0+j*step;

            double func = DetectorRect(gridX, gridY, deltaAngle*view, detector, detectors, views) *
                          CircleFunc(gridX, gridY, 1.0);

            H += func*std::exp(-2.0*M_PI*1i*(point[0]*gridX+point[1]*gridY)/2.0);
        }
    }
    H /= 4.0;
}


void CTPixelIntegal(double pixelSize, double integralStep,double deltaAngle,
                    int view, int detector, double gridX, double gridY,
                    int detectors, int views,
                    std::complex<double>& H)
{
    for(int i=0;i<std::floor(pixelSize/integralStep);i++)
    {
        for(int j=0;j<std::floor(pixelSize/integralStep);j++)
        {
            double point[2];
            point[0] = gridX+(i+1)*integralStep;
            point[1] = gridY+(j+1)*integralStep;

            double func = DetectorRect(point[0], point[1], deltaAngle*view, detector, detectors, views) *
                          PixelFunc(point[0], point[1], gridX, gridY, pixelSize) *
                          CircleFunc(point[0], point[1], 1.0);

            H += func*integralStep*integralStep;
        }
    }
    H /= pixelSize*pixelSize;
}


void CTDeltaIntegal(double pixelSize, double integralStep,double deltaAngle,
                    int view, int detector, double gridX, double gridY,
                    int detectors, int views,
                    std::complex<double>& H)
{
    double upperLeftX = gridX, upperLeftY = gridY;
    double upperRightX = gridX + pixelSize, upperRightY = gridY;
    double lowerLeftX = gridX, lowerLeftY = gridY + pixelSize;
    double lowerRightX = gridX + pixelSize, lowerRightY = gridY + pixelSize;

    double func1 = DetectorRect(upperLeftX, upperLeftY, deltaAngle*view, detector, detectors, views) *
                   CircleFunc(upperLeftX, upperLeftY, 1.0);

    double func2 = DetectorRect(upperRightX, upperRightY, deltaAngle*view, detector, detectors, views) *
                   CircleFunc(upperRightX, upperRightY, 1.0);

    double func3 = DetectorRect(lowerLeftX, lowerLeftY, deltaAngle*view, detector, detectors, views) *
                   CircleFunc(lowerLeftX, lowerLeftY, 1.0);

    double func4 = DetectorRect(lowerRightX, lowerRightY, deltaAngle*view, detector, detectors, views) *
                   CircleFunc(lowerRightX, lowerRightY, 1.0);

    if(func1 || func2 || func3 || func4)
        H = 1.0;
}


#endif // CTFUNC_H
