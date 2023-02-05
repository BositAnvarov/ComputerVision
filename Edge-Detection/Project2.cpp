#define _USE_MATH_DEFINES 
#include "Image.h"
#include <iostream>
#include <math.h>
#include <cmath>
using namespace std;

// Preconditions:  input image passed in is correctly allocated/formatted
// Postconditions: returns a new image that is a floating point image
Image fPoint(const Image& input)
{
    // result image with same dimensions
    Image res(input.getRows(), input.getCols());

    for (int i = 0; i < input.getRows(); i++)
    {
        for (int j = 0; j < input.getCols(); j++)
        {
            // set the specific pixel the float value 
            res.setFloat(i, j, input.getPixel(i, j).grey);
        }
    }

    // return the floating point image
    return res;
}

// Preconditions:  input image passed in is correctly allocated/formatted and is a floating image
// Postconditions: returns a new image that is a grey version of the floating point image
Image greyImage(const Image& input)
{
    // empty result image with same dimensions
    Image res(input.getRows(), input.getCols());

    for (int i = 0; i < input.getRows(); i++)
    {
        for (int j = 0; j < input.getCols(); j++)
        {
            // set the result image to grey for each pixel by using the static_cast<byte> method
            res.setGrey(i, j, static_cast<byte>(input.getFloat(i, j)));
        }
    }
    // return the grey image
    return res;
}

// Preconditions:  paramaters passed in are x or y, and the height or width accordingly
// Postconditions: returns a new integer that is the appopriate value according to the conditions of
//                 either being larger than boundary, less than zero or none
float checkCoordinate(const float& coordinate, int boundary)
{
    // if the coordinate is outside of the given boundary set it back to the boundary
    if (coordinate > (boundary - 1))
    {
        return float(boundary - 1);
    } // if the coordinate is below zero put it back to zero
    else if (coordinate < 0)
    {
        return 0.0f;
    }
    else
    {
        return coordinate;
    }
}

// Preconditions:  paramters(an appropriate image, either smoothing or gradient kernels, the origin of the kernel)
// Postconditions: returns a new image that is convloved with the given kernel
Image convolve(const Image& input, const Image& kernel, const int& cX, const int& cY)
{
    Image res(input.getRows(), input.getCols());

    // for each image row in result image
    for (int i = 0; i < input.getRows(); i++)
    {
        // for each image column in result image
        for (int j = 0; j < input.getCols(); j++)
        {
            // set running total to zero 
            float runningTotal = 0.0f;

            // for each kernel row
            for (int kI = 0; kI < kernel.getRows(); kI++)
            {
                // for each kernel column
                for(int kJ = 0; kJ < kernel.getCols(); kJ++)
                {

                    int x = int(checkCoordinate(float(j + (cX - kJ)), float(input.getCols())));
                    int y = int(checkCoordinate(float(i + (cY - kI)), float(input.getRows())));

                    // add result to the running total
                    runningTotal += input.getFloat(y, x) * kernel.getFloat(kI, kJ);
                }
            }

            // set the result image pixel to value of running total and reset running total for next pixel
            res.setFloat(i, j, runningTotal);
            runningTotal = 0.0f;
        }
    }

    // return result convolved image 
    return res;
}

// Preconditions: paramters(Magnitude image, x value of either r or p, y value of either r or p)
// Postconditions: using the alpha beta function compute the bilinear interpolation on the magnitude image 
//                 using x and y values
float interpolation(const Image& Gmag, const float& xVal, const float& yVal)
{
    // check the boundary for every float x and y passed
    float x = checkCoordinate(xVal , Gmag.getCols());
    float y = checkCoordinate(yVal, Gmag.getRows());

    // setup alpha and beta
    float alpha = y - floor(y);
    float beta = x - floor(x);

    // alpha beta equation given in class
    float res = ((1 - alpha) * (1 - beta) * Gmag.getFloat(floor(y), floor(x))) +
        (alpha * (1 - beta) * Gmag.getFloat(ceil(y), floor(x))) + ((1 - alpha) * beta * Gmag.getFloat(floor(y), ceil(x))) +
        (alpha * beta * Gmag.getFloat(ceil(y), ceil(x)));

    // return the result of bilinear interpollated image
    return res;

}

// Preconditions:  paramters(gradient in x direction, gradient in y direction)
// Postconditions: returns a new image that is a magnitude of gradientX and gradientY 
//                 using the equation Gmag = sqrt( (Gx * Gx) + (Gy * Gy) )
Image magnitudeImage(const Image& gX, const Image& gY)
{
    // empty result image with same dimensions of gradient
    Image res(gX.getRows(), gX.getCols());

    for (int i = 0; i < gY.getRows(); i++)
    {
        for (int j = 0; j < gX.getCols(); j++)
        {
            // square the gradient values 
            int gXsquared = gX.getFloat(i, j) * gX.getFloat(i, j);
            int gYsquared = gY.getFloat(i, j) * gY.getFloat(i, j);

            // magnitude equation commented above for each pixel
            res.setFloat(i, j, sqrt(gXsquared + gYsquared));
        }
    }

    // return the magnitude image
    return res;
}

// Preconditions:  parameter(a smoothed image). The purpose of the function is to detect edges using non maxima suppression and interpolation
// Postconditions: returns a new image that is a byte image that is grey scaled and dark with detected edges being white 
Image edgeDetectionFunc(const Image& floatingImage)
{
    Image res(floatingImage.getRows(), floatingImage.getCols());

    // set edge kernels for x and y direction
    Image gradientY(3, 1);
    Image gradientX(1, 3);

    // [-1, 0, 1] gradient kernel for x
    gradientX.setFloat(0, 0, -1);
    gradientX.setFloat(0, 1, 0);
    gradientX.setFloat(0, 2, 1);
    /*  _____
        |1/4|
        |1/2|
        |1/4|
        -----
    */
    gradientY.setFloat(0, 0, -1);
    gradientY.setFloat(1, 0, 0);
    gradientY.setFloat(2, 0, 1);

    Image Gx = convolve(floatingImage, gradientX, 1, 0);
    Image Gy = convolve(floatingImage, gradientY, 0, 1);
    Image Gmag = magnitudeImage(Gx, Gy);

    for (int i = 0; i < res.getRows(); i++)
    {
        for (int j = 0; j < res.getCols(); j++)
        {
            // check if magnitude is less at least 10 in the pixel
            if (!(10 <= Gmag.getFloat(i, j)))
            {
                res.setFloat(i, j, 0.0);
            }
            // do non maxima suppression 
            else
            {

                float gX = (Gx.getFloat(i, j)) / (Gmag.getFloat(i, j));
                float gY = (Gy.getFloat(i, j)) / (Gmag.getFloat(i, j));

                // r = q + g
                float rX = j + gX;
                float rY = i + gY;

                // p = q - g
                float pX = j - gX;
                float pY = i - gY;

                // Interpolate p and r 
                float r = interpolation(Gmag, rX, rY);
                float p = interpolation(Gmag, pX, pY);
                float g = Gmag.getFloat(i, j);

                // check edge: Gmag > gR and Gmag > gP 
                // if the conditions meet then set the edge brightness to 255(brightest)
                // if it does not then set the edge to 0.0(black)
                if (g > r && g > p) {
                    res.setFloat(i, j, 255.0);
                }
                else
                {
                    res.setFloat(i, j, 0.0);
                }
            }
        }
    }

    // return the smooth image with detected edges
    return res;
}

// Preconditions:  input image passed in is correctly allocated/formatted and is a floating image
// Postconditions: returns a new image that is a grey version of the floating point image
Image smoothFunc(const Image& input, const int& sigma)
{
    // convert image to a floating image
    Image floatingImage = fPoint(input);

    // create a smoothing kernel for X and Y
    Image smoothY(3, 1);
    Image smoothX(1, 3);

    // store floating image into a temporary image 
    Image res = floatingImage;

    // [1/4, 1/2, 1/4] smoothing kernel for x
    smoothX.setFloat(0, 0, .25);
    smoothX.setFloat(0, 1, .5);
    smoothX.setFloat(0, 2, .25);
    /*  _____
        |1/4|
        |1/2|
        |1/4|
        -----
    */
    smoothY.setFloat(0, 0, .25);
    smoothY.setFloat(1, 0, .5);
    smoothY.setFloat(2, 0, .25);

    // loop over sigma amount of times
    int num = 0;
    while (sigma > num)
    {
        // convolve both y and x kernel according to the origin of both
        res = convolve(res, smoothX, 1, 0);
        res = convolve(res, smoothY, 0, 1);

        num += 1;
    }

    // return the result image
    return res;
}

// Preconditions:  input image passed in is correctly allocated/formatted and is a floating image
// Postconditions: returns a new image that is a grey version of the floating point image 
int main(int argc, char* argv[])
{
    // run the sigma from the batchfile
    int sigma = 0;
    sscanf_s(argv[1], "%d", &sigma);


    Image input("test2.gif");
    
    // create an image and make it smooth using the smoothFunc, according to sigma
    // create an image and make it detect edges using the smooth image within the edgeDetectionFunc
    Image smoothImage = smoothFunc(input, sigma);
    Image edgeDetection = edgeDetectionFunc(smoothImage);

    // turn the images back into byte from floating point
    Image byteSmooth = greyImage(smoothImage);
    Image byteEdge = greyImage(edgeDetection);

    // write the images in grey
    byteSmooth.writeGreyImage("smooth.gif");
    byteEdge.writeGreyImage("edges.gif");

    return 0;
}
