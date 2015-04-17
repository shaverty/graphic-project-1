#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <math.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

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
      double		Z[3];
      double		color[3][3];

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
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].color[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].color[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].color[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

void line(double x1, double x2, double y, unsigned char * buffer, int width, int height, double  colors[][3])
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
    
    
    
    /*
     for (int j = leftIdx ; j <= rightIdx ; j++)
        {
            if (j < 0 || j >= width)
                continue;
            double proportion;
            if (rightPos != leftPos)
                proportion = (((double)j)-leftPos) / (rightPos - leftPos);
            else
                proportion = 1.0;
            unsigned char rgb[3];
            rgb[0] = (unsigned char) ceil441(255.0*(leftRGB[0] + proportion*(rightRGB[0]-leftRGB[0])));
            rgb[1] = (unsigned char) ceil441(255.0*(leftRGB[1] + proportion*(rightRGB[1]-leftRGB[1])));
            rgb[2] = (unsigned char) ceil441(255.0*(leftRGB[2] + proportion*(rightRGB[2]-leftRGB[2])));
            double z = leftZ + proportion*(rightZ-leftZ);
            Assign(j, i, rgb, z);
        }*/
        
        
    int leftldx = ceil441(x2);
    int rightldx = floor441(x1);
    
    for (int i = ceil441(x2); i <= floor441(x1); i++)
    {
	
      
      //cout << "TOP" << endl;
	buffer[3*((int)y*width+i)] = (unsigned char) ceil441(255.0*colors[0][0]);
	buffer[3*((int)y*width+i)+1] = (unsigned char) ceil441(255.0*colors[0][1]);
	buffer[3*((int)y*width+i)+2] = (unsigned char) ceil441(255.0*colors[0][2]);
    }
    
    leftldx = ceil441(x1);
    rightldx = floor441(x2);
    for (int i = ceil441(x1); i <= floor441(x2); i++)
    {
	//cout << " BOTTOM" << endl;
	buffer[3*((int)y*width+i)] = (unsigned char) ceil441(255.0*colors[1][0]);
	buffer[3*((int)y*width+i)+1] = (unsigned char) ceil441(255.0*colors[1][1]);
	buffer[3*((int)y*width+i)+2] = (unsigned char) ceil441(255.0*colors[1][2]);
    }
  
}

void draw(double x1, double y1, double x2, double y2, double x3, double y3, unsigned char *buffer, int width, int heigh, double colors[][3]){
  
  ///y1 == highest y
  ///y3 == lowest y
  double slope_left, slope_right, left, right, left_b, right_b;
  double height = y3 - y1;
  //cout << "in draw" << endl;
  if (height == 0)
  {
    cout<< "ERROR: on "<< times << endl;
    cout << "ceil of " << y3 << ": " << ceil441(y3) << " floor of " << y1 << ": " << floor441(y1) << " = " << height << endl;
    exit(-1);
  }
  
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
    for (int i = ceil441(y1); i <= floor441(y3); i++)
    {
      //cout << "entered first: " << left << " " << right << " " << i 
      left = slope_left == 0 ? left :(i-left_b)/slope_left;
      right = slope_right == 0? right : (i-right_b)/slope_right;
      line(right,left,i, buffer,width, heigh, colors);
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
 
   /*
    std::vector<Triangle> triangle;
   Triangle t;
   t.X[0] = 100;
   t.Y[0] = 100;
   t.X[1] = 150;
   t.Y[1] = 700;
   t.X[2] = 800;
   t.Y[2] = 800;
   t.color[0] = 255;
   t.color[1] = 255;
   t.color[2] = 255;
   triangle.push_back(t);
   */
   
   
   
   
   
   
   
  cout << triangle.size() << endl;
   
   
   for (int i = 0; i<triangle.size(); i++)
   {
      int order[3];
      
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
      
      
      
      
      ///flat bottom
      if (My1 == My2)
      {
	if (Mx1 < Mx2)
	   draw(Mx3,My3,Mx2,My2,Mx1,My2, buffer,width, height,triangle[i].color);
	else
	  draw(Mx3,My3,Mx1,My2,Mx2,My2, buffer,width, height,triangle[i].color);
      }
      ///flat top
      else if (My2 == My3)
      {
	if (Mx2< Mx3)
	  draw(Mx1,My1,Mx2,My2,Mx3,My2,buffer,width,height,triangle[i].color);
	else
	  draw(Mx1,My1,Mx3,My2,Mx2,My2,buffer,width, height,triangle[i].color);
      }
      ///needs splitting
      else if (My1 != My2 && My1 != My3 && My2 != My3)
      {
	double xNew = (Mx3-Mx1) * ((My2-My1)/(My3-My1)) + Mx1;
	
	if (xNew < Mx2)
	{
	  draw(Mx1, My1, xNew, My2, Mx2, My2,		buffer,width,height,triangle[i].color);
	  draw(Mx3, My3, Mx2, My2, xNew, My2,		buffer,width,height,triangle[i].color);
	}
	else
	{
	  draw(Mx1, My1, Mx2, My2, xNew, My2,		buffer,width,height,triangle[i].color);
	  draw(Mx3, My3, Mx2, My2, xNew, My2,		buffer,width,height,triangle[i].color);
	}
	times ++ ;
	
      }
      
	
   }	
   
   
    // buffer[0]= 255;
     
   
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   //cout << times << endl;
   WriteImage(image, "puddles1");
}
