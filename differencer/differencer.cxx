#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPNGReader.h>

using std::cerr;
using std::endl;

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    image->SetNumberOfScalarComponents(3);
    image->SetScalarType(VTK_UNSIGNED_CHAR);
    image->AllocateScalars();

    return image;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInput(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

vtkImageData *
ReadImage(const std::string &filename)
{
   vtkPNGReader *reader = vtkPNGReader::New();
   reader->SetFileName(filename.c_str());
   reader->GetOutput()->Update();
   if (reader->GetOutput() == NULL || reader->GetOutput()->GetNumberOfCells() == 0)
   {
       cerr << "Not able to read image " << filename << endl;
       exit(1);
   }
   return reader->GetOutput(); // memory leak
}

int main()
{
   std::string file1 = "frame750.png";
   std::string file2 = "/home/sean/Projects/project1f/build/Myframe750.png";
   vtkImageData *image1 = ReadImage(file1);
   int dims1[3];
   image1->GetDimensions(dims1);
   int npixels1 = dims1[0]*dims1[1];
   vtkImageData *image2 = ReadImage(file2);
   int dims2[3];
   image2->GetDimensions(dims2);
   int npixels2 = dims2[0]*dims2[1];
   if (npixels1 != npixels2)
   {
       cerr << "Image 1 has " << npixels1 << " pixels, while image 2 has " << npixels2 << " ... can't be a match!" << endl;
       exit(1);
   }
   cerr << "Both images have " << npixels1 << " pixels, good." << endl;

   unsigned char *buffer1 =
     (unsigned char *) image1->GetScalarPointer(0,0,0);
   unsigned char *buffer2 =
     (unsigned char *) image2->GetScalarPointer(0,0,0);

   bool foundDifference = false;
   int numDifferent = 0;

   // compare pixel by pixel.  Slowly overwrite the pixels in image1 to be the difference
   // map.
   for (int i = 0 ; i < npixels1 ; i++)
   {
        if (buffer1[3*i+0] != buffer2[3*i+0] ||
            buffer1[3*i+1] != buffer2[3*i+1] ||
            buffer1[3*i+2] != buffer2[3*i+2])
        {
          //cout << "pixel: " << i << ", original colors: " << (int)buffer1[3*i+0] << "/"<< (int)buffer1[3*i+1] <<"/"<< (int)buffer1[3*i+2] << ". My colors:  " << (int)buffer2[3*i+0] << "/"<<(int)buffer2[3*i+1] <<"/"<< (int)buffer2[3*i+2] << endl;
	  numDifferent++;
            buffer1[3*i+0] = 255;
            buffer1[3*i+1] = 255;
            buffer1[3*i+2] = 255;
	    //exit(0);
        }
        else
        {
            buffer1[3*i+0] = 0;
            buffer1[3*i+1] = 0;
            buffer1[3*i+2] = 0;
        }
   }

   cerr << "The number of different pixels is " << numDifferent << endl;
   if (numDifferent == 0)
      cerr << "Congratulations" << endl;
   else
   {
      cerr << "Writing file differenceMap.png" << endl;
      WriteImage(image1, "differenceMap");
   }
}
