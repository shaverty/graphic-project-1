#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <math.h>

using std::cerr;
using std::endl;
int times = 0;

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

void line(double x1, double x2, double y, unsigned char * buffer, int width, int height, unsigned char  colors[])
{
    if (x1<0)
      x1=0;
    if (x2 < 0)
      x2=0;
    if (y<0) 
      y=0;
    
    if (x1>=width) x1=width-1;
    if (x2>=width) x2=width-1;
    if (y>= height) y=height-1;
    //cout << "in line: " << x1<< " " << x2  << " "<< y << endl;
   
    
    for (int i = ceil441(x1); i <= floor441(x2); i++)
    {
	buffer[3*((int)y*width+i)] = colors[0];
	buffer[3*((int)y*width+i)+1] = colors[1];
	buffer[3*((int)y*width+i)+2] = colors[2];
    }
  
}

void draw(double x1, double y1, double x2, double y2, double x3, double y3, unsigned char *buffer, int width, int heigh, unsigned char colors[]){
  
  ///y1 == highest y
  ///y3 == lowest y
  double slope_left, slope_right, left, right, left_b, right_b;
  double height = y3 - y1;
  //cout << "in draw" << endl;
  if (height == 0)
  {
    cout<< "ERROR" << endl;
    exit(-1);
  }
  if (times == 99)
   cout << "poop" << endl;
  
  //slope_left = (x2-x1)/height; //inverse
  //slope_right = (x3-x1)/height; //inverse
  
 
  slope_left = (x2-x1) == 0 ? 0 : height/(x2-x1);
  slope_right = (x3-x1)== 0 ? 0 : height/(x3-x1);
  left_b = y1-x1*slope_left;
  right_b = y1-x1*slope_right;
  
  left = right = x1;
  //cout << y1 << " " << y3 << endl;
  if (y1<y3)
  {
    for (int i = ceil441(y1); i <= floor441(y3); ++i)
    {
      cout << "OMGGG" << endl;
      line(left, right, i, buffer,width,heigh, colors);
      left += ceil441(slope_left);
      right += floor441(slope_right);
    }
  }
  else
  {
    for (int i = floor441(y1); i>=ceil441(y3); i--)
    {
      left = slope_left == 0 ? left :(i-left_b)/slope_left;
      right = slope_right == 0? right : (i-right_b)/slope_right;
      line(right,left,i, buffer,width, heigh, colors);
      
    }
  }
  
}



int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;
   
   int width = 1000;
   int height = 1000;

   std::vector<Triangle> triangle = GetTriangles();
   
   
   
   
   
   
  
   
   
   for (int i = 0; i<triangle.size(); i++)
   {
      int order[3];
     //cout << "X coords: " << triangle[i].X[0] << " " << triangle[i].X[1]<< " " << triangle[i].X[2] <<endl;
   //cout << "Y coords: " << triangle[i].Y[0] << " " << triangle[i].Y[1]<< " " << triangle[i].Y[2] <<endl;
     
      
     //find min Y
     if (triangle[i].Y[0] < triangle[i].Y[1])
     {
	if (triangle[i].Y[0] < triangle[i].Y[2])
	  order[0] = 0;
	else
	   order[0] = 2;
     }
     else 
     {
       if (triangle[i].Y[1] < triangle[i].Y[2])
	   order[0] = 1;
       else
	  order[0] = 2;
     }
     
     //find max Y
     if (triangle[i].Y[0] > triangle[i].Y[1])
     {
       if (triangle[i].Y[0] > triangle[i].Y[2])
	 order[2] = 0;
       else
	 order[2] = 2;
     }
      else
      {
	if (triangle[i].Y[1] > triangle[i].Y[2])
	  order[2] = 1;
	else
	  order[2] = 2;
      }
      
      order[1] = 3-(order[0]+order[2]);
#define Mx1 triangle[i].X[order[0]]
#define Mx2 triangle[i].X[order[1]]
#define Mx3 triangle[i].X[order[2]]
#define My1 triangle[i].Y[order[0]]
#define My2 triangle[i].Y[order[1]]
#define My3 triangle[i].Y[order[2]]
     
      cout << My1 << " " << My2 << " " << My3 << endl;
      cout << Mx1 << " " << Mx2 << " " << Mx3 << endl;
      cout << endl;
      //cout << triangle[i].Y[order[0]] << " " << triangle[i].Y[order[2]]<< endl;
      if (My1 == My2)
      {
	//cout << "chose flat top" << endl;
	if (Mx1 < Mx2)
	   draw(Mx3,My3,Mx2,My2,Mx1,My2, buffer,width, height,triangle[i].color);
	else
	  draw(Mx3,My3,Mx1,My2,Mx2,My2, buffer,width, height,triangle[i].color);
      }
      else if (My2 == My3)
      {
	if (Mx2< Mx3)
	  draw(Mx1,My1,Mx2,My2,Mx3,My2,buffer,width,height,triangle[i].color);
	else
	  draw(Mx1,My1,Mx3,My2,Mx2,My2,buffer,width, height,triangle[i].color);
      }
	times ++;
   }	
   
   
    // buffer[0]= 255;
     
   
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   WriteImage(image, "allTriangles");
}
