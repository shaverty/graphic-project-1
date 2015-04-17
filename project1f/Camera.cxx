#include "Matrix.cxx"
class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3]; 	
    double          focus[3];
    double          up[3];		

    double 	    U[3];		//v1
    double	    V[3];		//v2
    double	    W[3];		//v3
    
    Matrix          ViewTransform();
    Matrix          CameraTransform();
    Matrix          DeviceTransform();
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

double dotProduct(double n[], double b[])
{
  double ret =  (n[0]*b[0]) + (n[1]*b[1]) + (n[2]*b[2]);
  return ret;
}

void normalize(double a[], double b[])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    
    b[0] = a[0]/length;
    b[1] = a[1]/length;
    b[2] = a[2]/length;
}

void crossProduct1D(double a[], double b[], double c[])
{
    double ret[3];
    ret[0] = a[1]*b[2]-a[2]*b[1];
    ret[1] = a[2]*b[0]-a[0]*b[2];
    ret[2] = a[0]*b[1]-a[1]*b[0];
    
    normalize(ret, ret);
    
    for (int i = 0; i<3; i++)
      c[i] = ret[i];
}

void makeFrame(Camera &c)
{
      crossProduct1D(c.up, c.position, c.U); 	//get v1 or U
      crossProduct1D(c.position, c.U, c.V);  	//get v2 or V
      normalize(c.position, c.W); 	 		//get v3 or W
}

Matrix
Camera::CameraTransform()
{
      Matrix m;
      double t[3];
      
      for (int i = 0 ; i < 3; i++){	//first column
	  m.A[i][0] =  U[i];
	  t[i] = 0.0-position[i];
      }
      m.A[3][0] = dotProduct(U,t); 
      
      
      for (int i = 0; i<3; i++)		//second column
	m.A[i][1] = V[i];
      m.A[3][1] = dotProduct(V,t);
      
      
      for (int i = 0; i<3; i++)		//third column
	m.A[i][2] = W[i];
      m.A[3][2] = dotProduct(W,t);
      
     for (int i = 0; i<4; i++)		//fourth column
	m.A[i][3] = i != 3? 0.0 : 1.0;
      
      
      return m;	
}

Matrix
Camera::ViewTransform()
{
      Matrix m;
      double cot = 1/tan(angle/2);
      
      for (int i = 0; i<4; i++)
	for (int j=0;j<4;j++)
	  m.A[i][j] = 0.0;
	
      m.A[0][0] = cot;
      m.A[1][1] = cot;
      m.A[2][2] = (far+near)/(far-near);
      m.A[2][3] = -1.0;
      m.A[3][2] = 2*(far*near)/(far-near);
      
      return m;
}

Matrix
Camera::DeviceTransform()
{
      Matrix m; 
      for (int i = 0; i<4; i++)
	for (int j=0;j<4;j++)
	  m.A[i][j] = 0.0;
	
      m.A[0][0] = 500.0;
      m.A[1][1] = 500.0;
      m.A[2][2] = 1.0;
      m.A[3][0] = 500.0;
      m.A[3][1] = 500.0;
      m.A[3][3] = 1.0;
      
      
      
      //return B.ComposeMatrices(m, B);
      return m;
}

