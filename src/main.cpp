#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fstream>


int SaveAsGadgetronRaw(std::string filename, float *raw, std::vector<int> dim)
{
    std::fstream file;

    file.open(filename, std::ios::out | std::ios::binary);

    if(!file.is_open())
        return -1;

    int ndim = dim.size();
    file.write((char*) &ndim, sizeof(int));
    long sizeOfRaw = 1;
    for(auto d : dim)
    {
        file.write((char*) &d, sizeof(int));
        sizeOfRaw *= d;
    }

    file.write((char*) raw, sizeOfRaw*sizeof(float));
}


auto rect = [](double value)
{
    if(value<=0.5 && value>=-0.5)
        return 1.0;
    else
        return 0.0;
};

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

auto PixelFunc = [](double x, double y, double pixelX, double pixelY, double pixelSize)
{
    if(x<pixelX+pixelSize && x>pixelX &&
       y<pixelY+pixelSize && y>pixelY)
        return 1.0;
    else
        return 0.0;
};

auto CircleFunc = [](double x, double y)
{
    std::complex<double> r(x, y);
    if(std::abs(r)<1)
        return 1.0;
    else
        return 0.0;
};


void PixelIntegal(double pixelSize, double integralStep,double deltaAngle,
                  int view, int detector, double gridX, double gridY,
                  double& H)
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
                          CircleFunc(point[0], point[1]);

            H += func*integralStep*integralStep;
        }
    }
    H /= pixelSize*pixelSize;
}



int main()
{
    int detectors = 32;
    double pixelSize = 2.0/128.0;
    double integralStep = pixelSize/50.0;
    int views = 16;
    double deltaAngle = (double)M_PI/views;

    double* gridX = new double[128];
    double* gridY = new double[128];

    for (int i=0;i<128;i++)
    {
        gridX[i] = -1.0 + pixelSize*i;
        gridY[i] = -1.0 + pixelSize*i;
    }

    double* H = new double[views*detectors*128*128];

    for(int i=0;i<views*detectors*128*128;i++)
    {
        H[i] = 0;
    }


    for(int a=0;a<views;a++)
    {
        for(int l=0;l<detectors;l++)
        {
            for(int x=0;x<128;x++)
            {
                for(int y=0;y<128;y++)
                {
                    PixelIntegal(pixelSize, integralStep, deltaAngle,
                                 a, l, gridX[x], gridY[y],
                                 H[y + x*128 + l*128*128 + a*detectors*128*128]);
                }
            }
            std::cout<<"l: "<<l<<std::endl;
        }
        std::cout<<"a: "<<a<<std::endl;
    }


    float* temp = new float[128*128*detectors*views];

    for(int i=0;i<128*128*detectors*views; i++)
    {
        temp[i] = H[i];
    }

    std::vector<int> dim =
    {
        128*128, detectors*views
    };
    SaveAsGadgetronRaw("H.real", temp, dim);


    delete [] gridX;
    delete [] gridY;
    delete [] H;

    return 0;
}
