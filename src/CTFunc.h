#include <complex>

#include "basicFunc.h"


auto DetectorRect = [&](double x, double y, double theta, int l)
{
    double rvec[2], nvec[2];
    rvec[0] = x;
    rvec[1] = y;
    nvec[0] = -std::sin(theta);
    nvec[1] = std::cos(theta);

    double dot = rvec[0]*nvec[0] + rvec[1]*nvec[1];

    return rect((dot-((l-16)/16.0 + 1.0/32.0))/(1.0/16.0));
};


void FourierIntegal(double pixelSize, double integralStep,double deltaAngle,
                    int view, int detector, double gridX, double gridY,
                    std::complex<double>& H)
{
    std::complex<double> i1 = std::complex<double>(0, 1);
    double step = 2.0/128.0;

    for(int i=0;i<128;i++)
    {
        for(int j=0;j<128;j++)
        {
            double point[2];
            point[0] = -1.0+i*step;
            point[1] = -1.0+j*step;

            double func = DetectorRect(gridX, gridY, deltaAngle*view, detector) *
                          CircleFunc(gridX, gridY, 1.0);

            H += func*std::exp(-2.0*M_PI*i1*(point[0]*gridX+point[1]*gridY)/2.0);
        }
    }
    H /= 4.0;
}


void PixelIntegal(double pixelSize, double integralStep,double deltaAngle,
                  int view, int detector, double gridX, double gridY,
                  std::complex<double>& H)
{
    for(int i=0;i<std::floor(pixelSize/integralStep);i++)
    {
        for(int j=0;j<std::floor(pixelSize/integralStep);j++)
        {
            double point[2];
            point[0] = gridX+(i+1)*integralStep;
            point[1] = gridY+(j+1)*integralStep;

            double func = DetectorRect(point[0], point[1], deltaAngle*view, detector) *
                          PixelFunc(point[0], point[1], gridX, gridY, pixelSize) *
                          CircleFunc(point[0], point[1], 1.0);

            H += func*integralStep*integralStep;
        }
    }
    H /= pixelSize*pixelSize;
}


void DeltaIntegal(double pixelSize, double integralStep,double deltaAngle,
                  int view, int detector, double gridX, double gridY,
                  std::complex<double>& H)
{
    double upperLeftX = gridX, upperLeftY = gridY;
    double upperRightX = gridX + pixelSize, upperRightY = gridY;
    double lowerLeftX = gridX, lowerLeftY = gridY + pixelSize;
    double lowerRightX = gridX + pixelSize, lowerRightY = gridY + pixelSize;

    double func1 = DetectorRect(upperLeftX, upperLeftY, deltaAngle*view, detector) *
                   CircleFunc(upperLeftX, upperLeftY, 1.0);

    double func2 = DetectorRect(upperRightX, upperRightY, deltaAngle*view, detector) *
                   CircleFunc(upperRightX, upperRightY, 1.0);

    double func3 = DetectorRect(lowerLeftX, lowerLeftY, deltaAngle*view, detector) *
                   CircleFunc(lowerLeftX, lowerLeftY, 1.0);

    double func4 = DetectorRect(lowerRightX, lowerRightY, deltaAngle*view, detector) *
                   CircleFunc(lowerRightX, lowerRightY, 1.0);

    if(func1 || func2 || func3 || func4)
        H == 1.0;
}
