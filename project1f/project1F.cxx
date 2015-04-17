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
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

#include <vtkDataSetWriter.h>
#include <stdio.h>
#include <stdlib.h>
#include "Camera.cxx"


using std::cerr;
using std::endl;
int times = 0;
using std::min;
using std::max;


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
      double 		normals[3][3];

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
    rdr->SetFileName("proj1e_geometry.vtk");
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
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
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
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0]; // no more multiple by 500...
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];


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

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
	 V[0]= 0;
	 V[1]= 1000;
	 V[2]= -1;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    }
  
    double V[3];	//viewing Direction;
    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;














int loopthing [4];
int NUMBER = 0;

Matrix transformTrianglesToDeviceSpace(){
  
      //Camera c = GetCamera(0,1000);
      //Camera c = GetCamera(250,1000);
      //Camera c = GetCamera(500,1000);
      Camera c = GetCamera(loopthing[NUMBER],1000);
      makeFrame(c);
       
      double views [3];
      
      normalize(c.position, views);
      for (int j=0; j<3; j++)
	lp.V[j] = -1.0*views[j];

      
      Matrix camera = c.CameraTransform();
      Matrix view = c.ViewTransform();
      Matrix device = c.DeviceTransform();
      
      view = view.ComposeMatrices(camera, view);
      device = view.ComposeMatrices(view, device);
      
      return device;
}

double CalculateShading(double normal[])
{
  double R[3], diffuse, specular, dot, shading; 
  dot = dotProduct(normal, lp.lightDir);
  
  for (int j =0; j< 3; j++)
    R[j] = 2*dot*normal[j]-lp.lightDir[j];
  
  diffuse = dot < 0? dot*-1.0*lp.Kd : dot*lp.Kd;
  
  specular = max(0.0,pow(dotProduct(R,lp.V),lp.alpha)*lp.Ks);
  
 
 shading =  lp.Ka + diffuse + specular;
  
  //cout << R[0] << "," << R[1]<< "," << R[2]<<"), with ambient = " << lp.Ka<< ", diffuse = "<< diffuse << ", specular = " << specular << ", for a total shading of " << shading<< endl;
  
  return shading;
}





void line(double x1, double x2, double y, Screen screen, double leftZ, double rightZ, double leftRGB[], double rightRGB[], double leftNorm[], double rightNorm[])
{
   if (y<0 || y >= screen.height)
     return;
        
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
      
	//double z = leftZ + ((i - x2) * zSlope);	
	double z = leftZ + proportion*(rightZ-leftZ);
	
	
 

	
	if (z > screen.z_buffer[(int)y*screen.width+i])
	{
      //cout << "TOP" << endl;
	  double norm[3], shading;
	  for (int j = 0; j < 3; j++){
	  norm[j] = leftNorm[j] + proportion*(rightNorm[j]-leftNorm[j]);
	  }
	  
	  shading = CalculateShading(norm);
      
      
	  unsigned char rgb[3];
	  for (int j = 0; j< 3; j++){
	    double thing = leftRGB[j] + proportion*(rightRGB[j]-leftRGB[j]);
	    thing = min(1.0, thing*shading);
	    rgb[j] = (unsigned char) ceil441(255.0*thing);
	  }
	
	  
	  screen.buffer[3*((int)y*screen.width+i)] = rgb[0];
	  screen.buffer[3*((int)y*screen.width+i)+1] = rgb[1];
	  screen.buffer[3*((int)y*screen.width+i)+2] = rgb[2];
	  
	  screen.z_buffer[(int)y*screen.width+i] = z;
  
	  //cout << "(" << i << ", " << y << ") with normal (" << norm[0] << ", " << norm[1] << ", " << norm[2] << ") and reflection = ("; 

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
       
       //double zSlope = (rightZ-leftZ) / (x2-x1);
       
      double proportion;
	if (x1 != x2)
	  proportion = (((double)i)-x2)/(x1-x2);
	else
	  proportion =1.0;
	
	
	//double z = leftZ + ((i - x1) * zSlope);
	double z = leftZ + proportion*(rightZ-leftZ);
	
	
	double norm[3], shading;
	for (int j = 0; j < 3; j++){
	 norm[j] = leftNorm[j] + proportion*(rightNorm[j]-leftNorm[j]);
	}
	
	shading = CalculateShading(norm);
	
	if (z > screen.z_buffer[(int)y*screen.width+i])
	{
	  //cout << "Bottom" << endl;
	  unsigned char rgb[3];
	  for (int j = 0; j< 3; j++){
	    double thing = leftRGB[j] + proportion*(rightRGB[j]-leftRGB[j]);
	    thing = min(1.0, thing*shading);
	    rgb[j] = (unsigned char)ceil441(255.0*thing);
	  }
	  
	  screen.buffer[3*((int)y*screen.width+i)] = rgb[0];
	  screen.buffer[3*((int)y*screen.width+i)+1] = rgb[1];
	  screen.buffer[3*((int)y*screen.width+i)+2] = rgb[2];
	  
	  screen.z_buffer[(int)y*screen.width+i] = z;
	  
	 // cout << "(" << i << ", " << y << ") with normal (" << norm[0] << ", " << norm[1] << ", " << norm[2] << ") and reflection = ("; 
	}
    }
   }
}



void draw(double x1, double y1, double z1, double color1[], double normal1[],
	  double x2, double y2, double z2, double color2[], double normal2[],
	  double x3, double y3, double z3, double color3[], double normal3[],
	  Screen screen){
  
  ///*1 == Point of triangle.  [max (flat bottom)]
  ///*3 == Base.		[min]
  ///Flat Top:    	x2 < x3
  ///Flat Bottom:	x3 < x2
  
  double slope_left, slope_right, left, right, left_b, right_b, za_slope, zb_slope, za, zb, leftRGB[3], rightRGB[3], leftNorm[3], rightNorm[3];
  double height = y3 - y1;
  
  if (height == 0)
  {
    cout<< "ERROR: on "<< times << endl;
    cout << "ceil of " << y3 << ": " << ceil441(y3) << " floor of " << y1 << ": " << floor441(y1) << " = " << height << endl;
    exit(-1);
  }
  

       
  
  zb_slope = (z3-z1) == 0? 0 : (z3-z1)/height; //slope for za
  za_slope = (z2-z1) == 0? 0 : (z2-z1)/height; //slope for zb
  
  
  slope_left = (x3-x1)== 0 ? 0 : height/(x3-x1);
  slope_right = (x2-x1) == 0 ? 0 : height/(x2-x1);
  left_b = y1-x1*slope_left;
  right_b = y1-x1*slope_right;
  
  left = right = x1;
  za = zb = z1;
  for (int i= 0; i < 3; i++){
    leftRGB[i] = rightRGB[i] = color1[i];
    leftNorm[i] = rightNorm[i] = normal1[i];
  }
  
  
  if (y1<y3){ //top flat triangle
  
    for (int i = ceil441(y1); i <= floor441(y3); i++)
    {
      left = slope_left == 0 ? left :(i-left_b)/slope_left;
      right = slope_right == 0? right : (i-right_b)/slope_right;
      
      za = za_slope == 0 ? za : z2 + ((i - y3) * za_slope);
      zb = zb_slope == 0 ? zb : z3 + ((i-y3) * zb_slope);

      
      
      for (int j = 0; j < 3; j++){
	leftRGB[j] = color2[j]==color1[j]   ? leftRGB[j] : color2[j] + (((double)i-y3) * ((color2[j]-color1[j])/height));
	leftNorm[j]=normal2[j]==normal1[j] ? leftNorm[j] : normal2[j] + (((double)i-y3) * ((normal2[j]-normal1[j])/height));
      }
       for (int j=0; j< 3; j++){
	 rightRGB[j] = color3[j] == color1[j]  ? rightRGB[j] : color3[j] + ((i-y3) * ((color3[j]-color1[j])/height));
	 rightNorm[j] = normal3[j]==normal1[j] ? rightNorm[j] : normal3[j] + ((i-y3) * ((normal3[j]-normal1[j])/height));
	 
       }
      
      line(left,right,i, screen, za, zb,leftRGB, rightRGB, leftNorm, rightNorm);
    }
  } 
else{ //bottom flat triangle
    for (int i = floor441(y1); i>=ceil441(y3); i--)
    {
      left = slope_left == 0 ? left :(i-left_b)/slope_left;
      right = slope_right == 0? right : (i-right_b)/slope_right;
      
      za = za_slope == 0 ? za : z1 + ((i-y1) * za_slope);
      zb = zb_slope == 0 ? zb : z1 + ((i-y1) * zb_slope);  
      
      for (int j = 0; j < 3; j++){
	leftRGB[j] = color2[j]==color1[j] 	? leftRGB[j]: color2[j] + ((i-y3) * ((color2[j]-color1[j])/height));
	leftNorm[j] = normal2[j]==normal1[j] 	? leftNorm[j]: normal2[j] + ((i-y3) * ((normal2[j]-normal1[j])/height));
	
      }
      
       for (int j=0; j< 3; j++){
	 rightRGB[j] = color3[j]==color1[j]	? rightRGB[j] : color3[j] + ((i-y3) * ((color3[j]-color1[j])/height));
	 rightNorm[j] = normal3[j]==normal1[j] ? rightNorm[j] :normal3[j]+ ((i-y3) * ((normal3[j]-normal1[j])/height));
      }
	 
      
     
	 /*
      cout << "Scanline " << i << " with left intercept at " << left << ", depth " << zb << ", and color " ;
      for (int j =0; j < 3; j++)
	cout << rightRGB[j] << "/";
      
       cout << ", right intercept at " << right << ", depth " << za << ", and color ";
      for (int j =0; j < 3; j++)
	cout << leftRGB[j] << "/";
      cout << endl;
      /** */
      
      line(right,left,i, screen ,zb, za, rightRGB, leftRGB, rightNorm, leftNorm);
      
    }
  
}

}



int main()
{
  loopthing[0] =0;
   loopthing[1]=250;
   loopthing[2]=500;
   loopthing[3]=750;
  for (NUMBER = 0; NUMBER<4; NUMBER++){
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
   
   
      Matrix matrix = transformTrianglesToDeviceSpace();
    
    
    
    //std::vector<Triangle> triangle = triangles;
      cout << triangle.size() << endl;

    for (int i = 0; i<triangle.size(); i++)
    {
      
      double *ptIn = new double[4];
      double *ptOut = new double[4];
      
      
      for (int j =0; j< 3; j++){
	ptIn[0] = triangle[i].X[j];
	ptIn[1] = triangle[i].Y[j];
	ptIn[2] = triangle[i].Z[j];
	ptIn[3] = 1.0;
	
	matrix.TransformPoint(ptIn, ptOut);
	
	triangle[i].X[j] = ptOut[0]/ptOut[3];
	triangle[i].Y[j] = ptOut[1]/ptOut[3];
	triangle[i].Z[j] = ptOut[2]/ptOut[3];
      }
	delete ptIn;
	delete ptOut;

      
	int *order = (int*)malloc(3*sizeof(int));
	
      //cout << "WORKING ON TRIANGLE " << i << endl;

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
  double *normal1= triangle[i].normals[order[0]];
  double *normal2= triangle[i].normals[order[1]];
  double *normal3= triangle[i].normals[order[2]];

	
    Screen screen;
    screen.buffer = buffer;
    screen.width = width;
    screen.height = height;
    screen.z_buffer = z_buffer;
	

	//cout << endl;
	///flat bottom
	if (My1 == My2)  //if bottom Y is equal to mid Y 
	{
	  if (Mx1 < Mx2) 				//if bottom X is left
	    draw(Mx3, My3, z3_g, color3, normal3, 
		Mx2, My2, z2_g, color2, normal2, 
		Mx1, My2, z1_g, color1, normal1,  screen);
	  else					//else switch
	    draw(Mx3, My3, z3_g, color3, normal3,
		Mx1, My2, z1_g, color1, normal1,
		Mx2, My2, z2_g, color2, normal2,  screen);
	}
	///flat top
	else if (My2 == My3) //if top Y is equal to mid Y
	{
	  if (Mx2< Mx3)				//if mid X is left
	    draw(Mx1,My1,z1_g,color1,normal1,
		Mx2,My2,z2_g,color2,normal2,
		Mx3,My2,z3_g,color3,normal3,  	screen);
	  else					//else switch
	    draw(Mx1,My1,z1_g,color1,normal1,
		Mx3,My2,z3_g,color3,normal3,
		Mx2,My2,z2_g,color2,normal2,  	screen);
	}
	///needs splitting
	else if (My1 != My2 && My1 != My3 && My2 != My3)
	{ 
	  double xNew = (Mx3-Mx1) * ((My2-My1)/(My3-My1)) + Mx1;
	  double zNew = z3_g + ((My2-My3)*((z1_g - z3_g) / (My1 - My3)));
	  
	  double *colorNew= new double[3];
	  for (int i = 0; i < 3; i++)
	    colorNew[i] = color3[i] + ((My2-My3)*((color1[i]-color3[i])/(My1-My3)));

	  double *normalNew = new double[3];
	  for (int i = 0; i < 3; i++)
	    normalNew[i] = normal3[i]+ ((My2-My3)*((normal1[i]-normal3[i])/(My1-My3)));
	  
	  
	  if (xNew < Mx2)
	  {
	    //top-flat
	    draw(Mx1, My1, z1_g, color1, normal1,
		xNew,My2, zNew,colorNew,normalNew,
		Mx2, My2, z2_g, color2, normal2,	screen);
	    draw(Mx3, My3, z3_g, color3, normal3,
		Mx2, My2, z2_g, color2, normal2,
		xNew,My2, zNew,colorNew,normalNew,	screen);
	    //bottom-flat
	  }
	  else
	  {
	    //top-flat
	    draw(Mx1, My1, z1_g, color1, normal1,
		Mx2, My2, z2_g, color2, normal2,
		xNew,My2, zNew,colorNew,normalNew,	screen);
	    draw(Mx3, My3, z3_g, color3, normal3,
		Mx2, My2, z2_g, color2, normal2,
		xNew,My2, zNew,colorNew,normalNew,	screen);
	    //bottom-flat
	  }
	  
	  free(colorNew);
	  free(normalNew);
	}
	
	
      //times ++ ;
      free(order);
    }	
      

      
      char Picture_buffer [12];// = "myframe";
      sprintf(Picture_buffer, "myframe%i", loopthing[NUMBER]);
    cout << "buffer: "<< Picture_buffer << endl;
    WriteImage(image, Picture_buffer);
   }
}
