#include <iostream>
#include <string>
#include <vector>
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmdata/dcfilefo.h"
#include <opencv2/opencv.hpp>
#include "ConfidenceMaps2DFacade.h"
using std::string;
using std::vector;
void disp_matrix2d(unsigned char * array, int width, int height);


int main() {
    char* volume_name = "../1.2.276.0.7230010.3.1.4.2355632107.136.1480927146.447";
    DicomImage * pDicomImg = new DicomImage(volume_name);
    std::cout << "read dicom file ...." << volume_name <<std::endl;
    int height = pDicomImg->getHeight();
    int width = pDicomImg->getWidth();
    int bits_num = pDicomImg->getDepth();
    int depth = pDicomImg->getFrameCount();
    std::cout << "pixel bitdepth---" << bits_num << std::endl;
    std::cout << "file size is  ---" << depth << " " << height << " " << width << std::endl;
    unsigned char *arr;
    arr = (unsigned char*)pDicomImg->getOutputData(8, 0, 0);
    //disp_matrix2d(arr, width, height);
    //pDicomImg->writeBMP("demo12.bmp", 8, 0);
    vector<double> img2d(height*width);
    unsigned char max = 0;
    unsigned char min = 255;
    for(int i = 0; i < width*height; i++)
    {
        if (*(arr+i) > max)
        {
            max = arr[i];
        }
        if (*(arr+i) < min)
        {
            min = arr[i];
        }
    }


    std::cout << "Min Pixel value--" << (int)min << "Max Pixel vale---" << (int)max << std::endl;
    std::cout << "convert to col major" << std::endl;
    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            img2d[i*height + j] = arr[j*width + i];
            //img2d[j*width + i] = arr[j*width + i];
        }
    }

    ConfidenceMaps2DFacade conf2d = ConfidenceMaps2DFacade();
    // conf2d.setImage(img2d, height, width, 1.5, true);
    // string solverType = "Vienna-CG-GPU";
    string solverType = "Eigen-LLT";
    conf2d.setSolver(solverType, 10000);
    conf2d.setImage(img2d, height, width, 1.5);

    std::vector<double> map = conf2d.computeMap();
    std::cout << "convert to row major" << std::endl;
    double val = 0;
    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            // unsigned char pixel_value = static_cast<unsigned char>(map[i*height + j]*255.0);
            val = map[i*height + j]*255;
            if(val < 0)
            {
                val = 0;
            }
            else if(val > 255)
            {
                val = 255;
            }
            unsigned char pixel_value = val;
            arr[j*width + i] = pixel_value;
        }
    }
    disp_matrix2d(arr, width, height);

    return 0;
}

void disp_matrix2d(unsigned char * array, int width, int height)
{
    // Create Mat object from array
    const cv::Mat img(cv::Size(width, height), CV_8U, array);
    // Display in window and wait for key press
    cv::namedWindow("foobar");
    cv::imshow("foobar", img);
    cv::waitKey(0);
}