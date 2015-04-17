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
      float *z_buffer;
      double *colors;

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

void line(double x1, double x2, double y, Screen screen, double leftZ, double rightZ, double leftRGB[], double rightRGB[])
{
   if (y<0 || y >= screen.height)
     return;
   
    
    

    
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
    
  
  if (x2 < x1){  
    for (int i = ceil441(x2); i <= floor441(x1); i++)
    {
      if (i < 0 || i >= screen.width)
	continue;
      
      double zSlope = (rightZ-leftZ) / (x1-x2);
      
      
	double proportion;
	if (x1 != x2)
	  proportion = (((double)i)-x2)/(x1-x2);
	else
	  proportion =1.0;
      
	double z = leftZ + ((i - x2) * zSlope);	
	//double z = leftZ + proportion*(rightZ-leftZ);
	
	if (z > screen.z_buffer[(int)y*screen.width+i])
	{
      //cout << "TOP" << endl;
	  unsigned char rgb[3];
	  rgb[0] = (unsigned char) ceil441(255.0*(leftRGB[0] + proportion*(rightRGB[0]-leftRGB[0])));
	  rgb[1] = (unsigned char) ceil441(255.0*(leftRGB[1] + proportion*(rightRGB[1]-leftRGB[1])));
	  rgb[2] = (unsigned char) ceil441(255.0*(leftRGB[2] + proportion*(rightRGB[2]-leftRGB[2])));
	  
	  screen.buffer[3*((int)y*screen.width+i)] = rgb[0];
	  screen.buffer[3*((int)y*screen.width+i)+1] = rgb[1];
	  screen.buffer[3*((int)y*screen.width+i)+2] = rgb[2];
	  
	  screen.z_buffer[(int)y*screen.width+i] = z;
	  //if ((int)y*screen.width+i > 999000 || (int)y*screen.width+i < 10 )
	  //cout << (int)y*screen.width+i << endl;

	  //cout<< "Assigning pixel (" << i << ", " << y << ") with color " << (int)rgb[0] << ", " << (int)rgb[1] << ", " << (int)rgb[2] << " and depth " << z << endl; 
	}
    }
}
    
    
   if (x1 < x2){
    leftldx = ceil441(x1);
    rightldx = floor441(x2);
    for (int i = ceil441(x1); i <= floor441(x2); i++)
    {
       if (i < 0 || i >= screen.width)
	continue;
       
       double zSlope = (rightZ-leftZ) / (x2-x1);
       
      double proportion;
	if (x1 != x2)
	  proportion = (((double)i)-x2)/(x1-x2);
	else
	  proportion =1.0;
	
	
	double z = leftZ + ((i - x1) * zSlope);
	//double z = leftZ + proportion*(rightZ-leftZ);
	
	if (z > screen.z_buffer[(int)y*screen.width+i])
	{
	  //cout << "Bottom" << endl;
	  unsigned char rgb[3];
	  rgb[0] = (unsigned char) ceil441(255.0*(leftRGB[0] + proportion*(rightRGB[0]-leftRGB[0])));
	  rgb[1] = (unsigned char) ceil441(255.0*(leftRGB[1] + proportion*(rightRGB[1]-leftRGB[1])));
	  rgb[2] = (unsigned char) ceil441(255.0*(leftRGB[2] + proportion*(rightRGB[2]-leftRGB[2])));
	  
	  screen.buffer[3*((int)y*screen.width+i)] = rgb[0];
	  screen.buffer[3*((int)y*screen.width+i)+1] = rgb[1];
	  screen.buffer[3*((int)y*screen.width+i)+2] = rgb[2];
	  
	  screen.z_buffer[(int)y*screen.width+i] = z;
	  
	  //if ((int)y*screen.width+i > 999900 || (int)y*screen.width+i < 10 )
	  //cout << (int)y*screen.width+i << endl;
	  
	  //cout << "Assigning pixel (" << i << ", " << y << ") with color " << (int)rgb[0] << ", " << (int)rgb[1] << ", " << (int)rgb[2] << " and depth " << z << endl; 
	}
    }
   }
  
}

void draw(double x1, double y1, double z1, double color1[],  double x2, double y2, double z2, double color2[], double x3, double y3, double z3, double color3[], Screen screen){
  
  ///y1 == highest y (flat bottom)
  ///y3 == lowest y
  ///x3 == left point
  ///x2 == right point
  
  double slope_left, slope_right, left, right, left_b, right_b, za_slope, zb_slope, za, zb, leftRGB[3], rightRGB[3];
  double height = y3 - y1;
  
  if (height == 0)
  {
    cout<< "ERROR: on "<< times << endl;
    cout << "ceil of " << y3 << ": " << ceil441(y3) << " floor of " << y1 << ": " << floor441(y1) << " = " << height << endl;
    exit(-1);
  }
  
  /*
  cout << "\n(" << x1 << ", " << y1 << ", " << z1 <<")/(" << color1[0] << ", " << color1[1] << ", " << color1[2] << "), "
       << "(" << x2 << ", " << y2 << ", " << z2 <<")/(" << color2[0] << ", " << color2[1] << ", " << color2[2] << "), "
       << "(" << x3 << ", " << y3 << ", " << z3 <<")/(" << color3[0] << ", " << color3[1] << ", " << color3[2] << ")" << endl; 
       cout << endl;
     /** */  
       
  
  zb_slope = (z3-z1) == 0? 0 : (z3-z1)/height; //slope for za
  za_slope = (z2-z1) == 0? 0 : (z2-z1)/height; //slope for zb
  
  
  slope_left = (x3-x1)== 0 ? 0 : height/(x3-x1);
  slope_right = (x2-x1) == 0 ? 0 : height/(x2-x1);
  left_b = y1-x1*slope_left;
  right_b = y1-x1*slope_right;
  
  left = right = x1;
  za = zb = z1;
  for (int i= 0; i < 3; i++)
    leftRGB[i] = rightRGB[i] = color1[i];
  
  
  if (y1<y3){
  
    for (int i = ceil441(y1); i <= floor441(y3); i++)
    {
      left = slope_left == 0 ? left :(i-left_b)/slope_left;
      right = slope_right == 0? right : (i-right_b)/slope_right;
      
      za = za_slope == 0 ? za : z2 + ((i - y3) * za_slope);
      zb = zb_slope == 0 ? zb : z3 + ((i-y3) * zb_slope);
      
      for (int j = 0; j < 3; j++)
	leftRGB[j] = (color2[j]-color1[j])/height == 0? leftRGB[j] : color2[j] + ((i-y3) * ((color2[j]-color1[j])/height));
      
      /*
      cout << "Scanline " << i << " with left intercept at " << right << ", depth " << za << ", and color " ;
      for (int j =0; j < 3; j++)
	cout << leftRGB[j] << "/";
      /** */
      
      for (int j=0; j< 3; j++)
	 rightRGB[j] = (color3[j]-color1[j])/height == 0? leftRGB[j] : color3[j] + ((i-y3) * ((color3[j]-color1[j])/height));
      
      /*
      cout << ", right intercept at " << left << ", depth " << zb << ", and color ";
      for (int j =0; j < 3; j++)
	cout << rightRGB[j] << "/";  
      cout << endl;
      /** */
     /* 
     if (left > right)
     {
      double temp =0;
      temp = left;
      left = right;
      right = temp;
     }*/
      
      line(left,right,i, screen, za, zb,leftRGB, rightRGB);
    }
  }
  else
  {
    for (int i = floor441(y1); i>=ceil441(y3); i--)
    {
      left = slope_left == 0 ? left :(i-left_b)/slope_left;
      right = slope_right == 0? right : (i-right_b)/slope_right;
      
      za = za_slope == 0 ? za : z1 + ((i-y1) * za_slope);
      zb = zb_slope == 0 ? zb : z1 + ((i-y1) * zb_slope);
      
      for (int j = 0; j < 3; j++)
	leftRGB[j] = (color2[j]-color1[j])/height == 0? color1[j]: color2[j] + ((i-y3) * ((color2[j]-color1[j])/height));
      
         for (int j=0; j< 3; j++)
	 rightRGB[j] = (color3[j]-color1[j])/height == 0? leftRGB[j] : color3[j] + ((i-y3) * ((color3[j]-color1[j])/height));
      
     /*
      cout << "Scanline " << i << " with left intercept at " << left << ", depth " << zb << ", and color " ;
      for (int j =0; j < 3; j++)
	cout << rightRGB[j] << "/";
      /** */
      
   
       /*
       cout << ", right intercept at " << right << ", depth " << za << ", and color ";
      for (int j =0; j < 3; j++)
	cout << leftRGB[j] << "/";
      cout << endl;
      /** */
      /*
      if (left > right)
     {
      double temp =0;
      temp = left;
      left = right;
      right = temp;
     }
     */
      
      line(right,left,i, screen ,zb, za, rightRGB, leftRGB);
      
    }
  }
  

}



int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   float *z_buffer = (float*)malloc(npixels*sizeof(float));
   
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   for (int i = 0; i < npixels; i++)
      z_buffer[i]=-1.0;
   
   
   int width = 1000;
   int height = 1000;


   std::vector<Triangle> triangle = GetTriangles();
 
   
   // std::vector<Triangle> triangle; Triangle t;
   /*
   t.X[0] = 100;
   t.Y[0] = 100;
   t.X[1] = 150;
   t.Y[1] = 700;
   t.X[2] = 800;
   t.Y[2] = 800;
   
   (551.02, 877.551, -0.953298)/(0.117997, 0.117997, 0.569589), 
   (551.02, 864.703, -0.938775)/(0.100532, 0.100532, 0.538104), 
   (567.985, 864.703, -0.941226)/(0.103373, 0.103373, 0.543226)
   */
  
   
   
   
   
   
   
   
  cout << triangle.size() << endl;

   for (int i = 0; i<triangle.size(); i++)
   {
     //if (i == 19)
       //continue;
      int *order = (int*)malloc(3*sizeof(int));
      
     //cout << "WORKING ON TRIANGLE " << i << endl;
      //cout << endl ;
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
double Mx1 = triangle[i].X[order[0]];
double Mx2 = triangle[i].X[order[1]];
double Mx3 = triangle[i].X[order[2]];
double My1 = triangle[i].Y[order[0]];
double My2 = triangle[i].Y[order[1]];
double My3 = triangle[i].Y[order[2]];
double z1_g  = triangle[i].Z[order[0]];
double z2_g  = triangle[i].Z[order[1]];
double z3_g  = triangle[i].Z[order[2]];
double *color1 = triangle[i].color[order[0]];
double *color2 = triangle[i].color[order[1]];
double *color3 = triangle[i].color[order[2]];

      
   Screen screen;
   screen.buffer = buffer;
   screen.width = width;
   screen.height = height;
   screen.z_buffer = z_buffer;
      
      
      ///flat bottom
      if (My1 == My2)
      {
	if (Mx1 < Mx2)
	  draw(Mx3,My3,z3_g,color3, Mx2,My2,z2_g,color2, Mx1,My2,z1_g,color1,  screen);
	else
	  draw(Mx3,My3,z3_g,color3, Mx1,My2,z1_g,color1, Mx2,My2,z2_g,color2,  screen);
      }
      ///flat top
      else if (My2 == My3)
      {
	if (Mx2< Mx3)
	  draw(Mx1,My1,z1_g,color1, Mx2,My2,z2_g,color2, Mx3,My2,z3_g,color3,  screen);
	else
	  draw(Mx1,My1,z1_g,color1, Mx3,My2,z3_g,color3, Mx2,My2,z2_g,color2,  screen);
      }
      ///needs splitting
      else if (My1 != My2 && My1 != My3 && My2 != My3)
      {
	//double *xNew = new double; 
	 double xNew = (Mx3-Mx1) * ((My2-My1)/(My3-My1)) + Mx1;
	//double *zNew = new double;
	double zNew = z3_g + ((My2-My3)*((z1_g - z3_g) / (My1 - My3)));
	double *colorNew= new double[3];
	colorNew[0] = color3[0] + ((My2-My3)*((color1[0]-color3[0])/(My1-My3)));
	colorNew[1] = color3[1] + ((My2-My3)*((color1[1]-color3[1])/(My1-My3)));
	colorNew[2] = color3[2] + ((My2-My3)*((color1[2]-color3[2])/(My1-My3)));
	
	//cout << "\ncolor info: \n" << color1[0] << " " << color1[1] << " " << color1[2] << "\n"<< color2[0] << " " << color2[1] << " " << color2[2] << "\n" << color3[0] << " " << color3[1] << " " << color3[2] << endl;
	
	//cout << colorNew[0] << " " << colorNew[1] << " " << colorNew[2] << "\n" << endl;
	
	if (xNew < Mx2)
	{
	  //top-flat
	  draw(Mx1,My1,z1_g,color1,  xNew,My2,zNew,colorNew,  Mx2,My2,z2_g,color2,		screen);
	  draw(Mx3,My3,z3_g,color3,  Mx2,My2,z2_g,color2,  	 xNew,My2,zNew,colorNew,	screen);
	  //bottom-flat
	}
	else
	{
	  //top-flat
	  draw(Mx1,My1,z1_g,color1,  Mx2,My2,z2_g,color2, xNew,My2,zNew,colorNew,		screen);
	  draw(Mx3,My3,z3_g,color3,  Mx2,My2,z2_g,color2, xNew,My2,zNew,colorNew,		screen);
	  //bottom-flat
	}
	times ++ ;
	free(colorNew);
      }
      

 
    free(order);
   }	
   
   
    // buffer[0]= 255;
     
   
   


   //cout << times << endl;
   WriteImage(image, "puddles");
}
