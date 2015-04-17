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
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1c_geometry.vtk");
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
    vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
    float *color_ptr = colors->GetPointer(0);
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
        tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
        tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
        tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
    }

    return tris;
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
    
    
    for (int i = ceil441(x2); i <= floor441(x1); i++)
    {
	
      
      //cout << "TOP" << endl;
	buffer[3*((int)y*width+i)] = colors[0];
	buffer[3*((int)y*width+i)+1] = colors[1];
	buffer[3*((int)y*width+i)+2] = colors[2];
    }
    
    
    for (int i = ceil441(x1); i <= floor441(x2); i++)
    {
	//cout << " BOTTOM" << endl;
	buffer[3*((int)y*width+i)] = colors[0];
	buffer[3*((int)y*width+i)+1] = colors[1];
	buffer[3*((int)y*width+i)+2] = colors[2];
    }
  
}

void draw(double x1, double y1, double x2, double y2, double x3, double y3, unsigned char *buffer, int width, int heigh, unsigned char colors[]){
  
  ///y1 == highest y (flat-bottom)
  ///y3 == lowest y 
  
  ///x3 == left point
  ///x2 == right point
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
  
  
  slope_left = (x3-x1) == 0 ? 0 : height/(x3-x1);
  slope_right = (x2-x1)== 0 ? 0 : height/(x2-x1);
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
      line(left,right,i, buffer,width, heigh, colors);
    }
  }
  else
  {
    for (int i = floor441(y1); i>=ceil441(y3); i--)
    {
      left = slope_left == 0 ? left :(i-left_b)/slope_left;
      right = slope_right == 0? right : (i-right_b)/slope_right;
      line(left,right,i, buffer,width, heigh, colors);
      
    }
  }
  

}



int main()
{
   vtkImageData *image = NewImage(1786, 1344);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1786*1344;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;
   
   int width = 1786;
   int height = 1344;


   //std::vector<Triangle> triangle = GetTriangles();
 
   
    std::vector<Triangle> triangle;
   Triangle t;
   
   t.X[0] = 100;
   t.Y[0] = 100;
   t.X[1] = 150;
   t.Y[1] = 700;
   t.X[2] = 800;
   t.Y[2] = 800;
   /*
   (551.02, 877.551, -0.953298)/(0.117997, 0.117997, 0.569589), (551.02, 864.703, -0.938775)/(0.100532, 0.100532, 0.538104), (567.985, 864.703, -0.941226)/(0.103373, 0.103373, 0.543226)

   */
   /*
   t.X[0] = 551.02;
   t.Y[0] = 877.551;  //-0.953298
   t.X[1] = 551.02;	
   t.Y[1] = 864.703;//0.569589
   t.X[2] = 567.985;
   t.Y[2] = 864.703; //0.941226
   */
   t.color[0] = 255;
   t.color[1] = 255;
   t.color[2] = 255;
   triangle.push_back(t);
   
   
   
   
   
   
   
   
  cout << triangle.size() << endl;
   
   
   for (int i = 0; i<triangle.size(); i++)
   {
      int order[3];
     /*
      cout << "X coords: " << triangle[i].X[0] << " " << triangle[i].X[1]<< " " << triangle[i].X[2] <<endl;
   cout << "Y coords: " << triangle[i].Y[0] << " " << triangle[i].Y[1]<< " " << triangle[i].Y[2] <<endl;
      cout << endl;
      triangle.size()
     */
      
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
     
      
      
      //cout << triangle[i].Y[order[0]] << " " << triangle[i].Y[order[2]]<< endl;
      
      
      
      
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
	
	cout << xNew << " " << Mx2 << endl;
	if (xNew < Mx2)
	{
	  ///draws top-flat triangle 			(bottom of split)
	  draw(Mx1, My1, xNew, My2, Mx2, My2,		buffer,width,height,triangle[i].color);
	  cout << Mx1 << " " << My1 << " " <<xNew<< " " << My2<< " " <<Mx2<< " " <<My2 << endl;
	  ///draws bottom-flat triangle 		(top of split)
	  draw(Mx3, My3, Mx2, My2, xNew, My2,		buffer,width,height,triangle[i].color);
	  cout << Mx3 << " " << My3 << " " <<Mx2<< " " << My2<< " " <<xNew<< " " <<My2 << endl;
	}
	else
	{
	  ///draws top-flat triangle 			(bottom of split)
	  draw(Mx1, My1, Mx2, My2, xNew, My2,		buffer,width,height,triangle[i].color);   
	  
	  ///draws bottom-flat triangle 		(top of split)
	  draw(Mx3, My3, Mx2, My2, xNew, My2,		buffer,width,height,triangle[i].color);
	}
	times ++ ;
	
      }
      
	
   }	
   
   
    // buffer[0]= 255;
     
   
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1786;
   screen.height = 1344;

   //cout << times << endl;
   WriteImage(image, "puddles");
}
