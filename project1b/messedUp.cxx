#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <math.h>

using std::cerr;
using std::endl;

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}


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

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      unsigned char color[3];

  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;

  // would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
   std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       if (i == 50)
           rv[i].X[firstPt] = -10;
       rv[i].Y[firstPt] = posJ;
       rv[i].X[(firstPt+1)%3] = posI+99;
       rv[i].Y[(firstPt+1)%3] = posJ;
       rv[i].X[(firstPt+2)%3] = posI+i;
       rv[i].Y[(firstPt+2)%3] = posJ+10*(idxJ+1);
       if (i == 95)
          rv[i].Y[(firstPt+2)%3] = 1050;
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
   }

   return rv;
}

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangle = GetTriangles();
   
   
   for (int i = 0; i<99; i++)
   {
     float invslope1 = abs((triangle[i].X[0] - triangle[i].X[2]) / (triangle[i].Y[0]-triangle[i].Y[2]));
     float invslope2 = abs((triangle[i].X[1] - triangle[i].X[2]) / (triangle[i].Y[1]-triangle[i].Y[2]));
     
     float curx1 = triangle[i].X[2];
     float curx2 = triangle[i].X[2];
     for (int scanY = triangle[i].Y[2]; scanY >= triangle[i].Y[0]; scanY--)
     {
	for (int j = (int)ceil(curx1); j < (int)floor(curx2); j++)
	{
	  if (scanY == 0)
	  {
	   buffer[j*3] = triangle[i].color[0];
	   buffer[j*3 +1] = triangle[i].color[1];
	   buffer[j*3 +2] = triangle[i].color[2];
	  }
	  else
	  {
	    buffer[3*(scanY*1000+j)] = triangle[i].color[0];
	    buffer[3*(scanY*1000+j)+1] = triangle[i].color[1];
	    buffer[3*(scanY*1000+j)+2] = triangle[i].color[2];
	    
	    //cout << "things: " << 3*(scanY*1000+j) << " "<< 3*(scanY*1000+j)+1 << " " << 3*(scanY*1000+j)+2 << endl;
	  }
	}
	curx1 +=invslope1;
	curx2 +=invslope2;
	//cout << curx1 << " " << curx2 << endl;
	//cout << scanY << endl;
     }
     
   }
   
   //cout << "X coords: " << triangle[0].X[0] << " " << triangle[0].X[1]<< " " << triangle[0].X[2] <<endl;
   //cout << "Y coords: " << triangle[0].Y[0] << " " << triangle[0].Y[1]<< " " << triangle[0].Y[2] <<endl;
     
   
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   WriteImage(image, "allTriangles");
}
