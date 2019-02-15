#ifndef BASICFUNC_H
#define BASICFUNC_H

#include <fstream>
#include <cmath>
#include <complex>
#include <vector>


auto rect = [](double value)
{
    if(value<=0.5 && value>=-0.5)
        return 1.0;
    else
        return 0.0;
};


auto PixelFunc = [](double x, double y, double pixelX, double pixelY, double pixelSize)
{
    if(x<pixelX+pixelSize && x>pixelX &&
       y<pixelY+pixelSize && y>pixelY)
        return 1.0;
    else
        return 0.0;
};

auto CircleFunc = [](double x, double y, double R)
{
    std::complex<double> r(x, y);
    if(std::abs(r)<R)
        return 1.0;
    else
        return 0.0;
};


template<class T> int SaveAsGadgetronRaw(std::string filename, T *raw, std::vector<int> dim)
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

    file.write((char*) raw, sizeOfRaw*sizeof(T));
    
    file.close();
    
    return 0;
}


#endif // BASICFUNC_H
