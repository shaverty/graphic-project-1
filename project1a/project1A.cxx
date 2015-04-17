#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

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


int main()
{
   std::cerr << "In main!" << endl;
   vtkImageData *image = NewImage(1024, 1350);
   image->Print(cerr);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   cerr << "buffer is " << (void *) buffer << endl;
   
   int curr = 0;
   for (int x=0; x < 27; x++)
   {
      int R = 0;
      int G = 0;
      int B = 0;
      
       ///red
    if (x/9 == 0)
      R=0;
    else if (x/9 == 1)
      R = 128;
    else if (x/9 ==2)
      R=255;
    ///red
    
    ///green
      if ((x/3) %3 ==0)
	G=0;
      else if ((x/3) %3 ==1)
	G=128;
      else if ((x/3) %3 ==2)
	G=255;
    ///green
    
    ///blue
      if (x%3 == 0)
	  B=0;
      else if (x%3 == 1)
	  B=128;
      else if (x%3 == 2)
	  B=255;
    ///blue
	    
	  
      for (int i=curr; i < curr+50*1024*3; i+=3){ 
	buffer[i]   = R;
	buffer[i+1] = G;
	buffer[i+2] = B;
      }
      
      curr += 50*1024*3;
   }
   
   WriteImage(image, "oneRedPixel");
}
