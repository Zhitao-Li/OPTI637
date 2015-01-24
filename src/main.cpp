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


int main()
{
    int detectors = 32;
    double pixelSize = 2.0/256.0;
    double integralStep = pixelSize/10.0;
    int views = 16;
    double deltaAngle = (double)M_PI/views;

    double* gridX = new double[256];
    double* gridY = new double[256];

    for (int i=0;i<256;i++)
    {
        gridX[i] = -1.0 + pixelSize*i;
        gridY[i] = -1.0 + pixelSize*i;
    }

    double CenterX[10] =
    {
        0, 0, 0.22, -0.22, 0, 0, 0, -0.08, 0, 0.06
    };
    double CenterY[10] =
    {
        0, -0.0184, 0, 0, 0.35, 0.1, -0.1, -0.605, -0.605, -0.605
    };
    double MajorAxis[10] =
    {
        0.92, 0.874, 0.31, 0.41, 0.25, 0.046, 0.046, 0.046, 0.023, 0.046
    };
    double MinorAxis[10] =
    {
        0.69, 0.6624, 0.11, 0.16, 0.21, 0.046, 0.046, 0.023, 0.023, 0.023
    };
    double RotationAngle[10] =
    {
        90, 90, 72, 108, 90, 0, 0, 0, 0, 90
    };
    for(int i=0;i<10;i++)
    {
        RotationAngle[i] = RotationAngle[i]/180.0*M_PI;
    }
    double FuncOffset[10]
    {
        2, -0.98, -0.2, -0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1
    };

    double* area = new double[views*detectors*256*256];
    double* accum = new double[views*detectors];

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


    for(int i=0;i<views*detectors*256*256;i++)
    {
        area[i] = 0;
    }


    for(int a=0;a<views;a++)
    {
        for(int l=0;l<detectors;l++)
        {
            for(int x=0;x<256;x++)
            {
                for(int y=0;y<256;y++)
                {
                    for(int i=0;i<std::floor(pixelSize/integralStep);i++)
                    {
                        for(int j=0;j<std::floor(pixelSize/integralStep);j++)
                        {
                            double point[2];
                            point[0] = gridX[x]+(i+1)*integralStep;
                            point[1] = gridY[y]+(j+1)*integralStep;

                            double func = DetectorRect(point[0], point[1], deltaAngle*a, l) *
                                          PixelFunc(point[0], point[1], gridX[x], gridY[y], pixelSize) *
                                          CircleFunc(point[0], point[1]);

                            if(func!=0)
                            {
                                for(int c=0;c<10;c++)
                                {
                                    double rotate[4];
                                    rotate[0] = cos(RotationAngle[c]);
                                    rotate[1] = -sin(RotationAngle[c]);
                                    rotate[2] = sin(RotationAngle[c]);
                                    rotate[3] = cos(RotationAngle[c]);

                                    double p[2] =
                                    {
                                        rotate[0]*point[0] + rotate[1]*point[1],
                                        rotate[2]*point[0] + rotate[3]*point[1]
                                    };

                                    if((p[0]-CenterX[c])*(p[0]-CenterX[c])/(MajorAxis[c]*MajorAxis[c]) +
                                       (p[1]-CenterY[c])*(p[1]-CenterY[c])/(MinorAxis[c]*MinorAxis[c]) <= 1.0)
                                        accum[l + a*detectors] += FuncOffset[c]*integralStep*integralStep;
                                }
                            }

                            area[y + x*256 + l*256*256 + a*detectors*256*256] += func*integralStep*integralStep;
                        }
                    }
                }
            }
            std::cout<<"l: "<<l<<std::endl;
        }
        std::cout<<"a: "<<a<<std::endl;
    }


    float* temp = new float[256*256*detectors*views];

    for(int i=0;i<256*256*detectors*views; i++)
    {
        temp[i] = area[i];
    }

    float* temp1 = new float[detectors*views];

    for(int i=0;i<detectors*views; i++)
    {
        temp1[i] = accum[i];
    }

    std::vector<int> dim =
    {
        256, 256, detectors, views
    };
    SaveAsGadgetronRaw("area.real", temp, dim);

    std::vector<int> dim1 =
    {
        detectors, views
    };
    SaveAsGadgetronRaw("accum.real", temp1, dim1);


    delete [] gridX;
    delete [] gridY;
    delete [] area;
    delete [] accum;

    return 0;
}
