#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <iostream>
#include <string>
#include <math.h>
#include <sys/stat.h>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;
#define NORMALS

double ceil(double f)
{
    return ceil(f-0.00001);
}

double floor(double f)
{
    return floor(f+0.00001);
}

vtkImageData* NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    return img;
}

void WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

double getEndValue(double x1, double x2, double y1, double y2, double y){

    if (x1 == x2){
        return x1;
    }
    else{
        //y = slope * x + n
        double slope = (y1-y2) / (x1-x2);
        double n = y1 - slope * x1;
        double end = (y-n) / slope;
        return end;
    }
}

void normalize(double* v, double* result, int size){
    double m = 0;
    for (int i=0; i<size; i++){
        m = m + v[i] * v[i];
    }
    m = sqrt(m);
    for (int j=0; j<size; j++){
        result[j] = v[j] / m;
    }
}

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

class Vertex
{
public:
    double X;
    double Y;
    double Z;
    double color[3];
    double normal[3];
    double shading;
    double AMB;
    double DIFF;
    double SPEC;

    void printVertex(){
        cout<<std::setprecision(10);
        cout<<X<<", "<<Y<<", "<<Z<<"\n";
        cout<<"color: ";
        cout<<ceil(color[0]*255)<<", ";
        cout<<ceil(color[1] * 255)<<", "<<ceil(color[2] * 255)<<"\n";
        cout<<"normals: ";
        cout<<normal[0]<<", "<<normal[1]<<", "<<normal[2]<<"\n";
        cout<<"Shading: "<<shading<<", "<<"AMB = "<<AMB<<", DIFF = "<<DIFF<<", SPEC = "<<SPEC<<"\n\n";
    }

    double getX(){
        return X;
    }

    double getY(){
        return Y;
    }

    double getZ(){
        return Z;
    }

    double* getColor(){
        return color;
    }

    double* getNormal(){
        return normal;
    }
    
    double getShading(){
        return shading;
    }

    void setXYZ(double newx, double newy, double newz){
        X = newx;
        Y = newy;
        Z = newz;
    }

    void setColor(double r, double g, double b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void setNormal(double n1, double n2, double n3){
        normal[0] = n1;
        normal[1] = n2;
        normal[2] = n3;
    }

    void setShading(double s){
        shading = s;
    }

    void CalculatePhongShading(LightingParameters &lp, double *viewDirection){

        AMB = lp.Ka;

        double LN = lp.lightDir[0] * normal[0] + lp.lightDir[1] * normal[1] + lp.lightDir[2] * normal[2];
        DIFF = lp.Kd * abs(LN);

        double R[3];
        for (int k=0; k<3; k++){
            R[k] = 2 * normal[k] * LN - lp.lightDir[k];
        }
        normalize(R, R, 3);
        double RV = R[0] * viewDirection[0] + R[1] * viewDirection[1] + R[2] * viewDirection[2];

        
        SPEC = fmax(0, lp.Ks * pow(RV, lp.alpha));
        
        shading = AMB + DIFF + SPEC;
    }

};

Vertex interpolate(Vertex v1, Vertex v2, Vertex v3){

    double t;

    if (v2.getY() != v1.getY()){
        //cout<<"if x\n";
        t = (v3.getY() - v1.getY()) / (v2.getY() - v1.getY());
    }
    else if (v2.getX() != v1.getX()){
        //cout<<"if y\n";
        t = (v3.getX() - v1.getX()) / (v2.getX() - v1.getX());
    }
    else{
        //cout<<"if none\n";
        t = 0;
    }

    double r = v1.getColor()[0] + t * (v2.getColor()[0] - v1.getColor()[0]);
    double g = v1.getColor()[1] + t * (v2.getColor()[1] - v1.getColor()[1]);
    double b = v1.getColor()[2] + t * (v2.getColor()[2] - v1.getColor()[2]);

    double n1 = v1.getNormal()[0] + t * (v2.getNormal()[0] - v1.getNormal()[0]);
    double n2 = v1.getNormal()[1] + t * (v2.getNormal()[1] - v1.getNormal()[1]);
    double n3 = v1.getNormal()[2] + t * (v2.getNormal()[2] - v1.getNormal()[2]);

    double s = v1.getShading() + t * (v2.getShading() - v1.getShading());

    v3.setColor(r, g, b);
    v3.setNormal(n1, n2, n3);
    v3.setShading(s);
    
    return v3;
}

class Triangle
{
  public:
    double         X[3];
    double         Y[3];
    double         Z[3];
    double         colors[3][3];
    double         normals[3][3];
    double         shadings[3];
    double         cPosition[3];
    Vertex         vertices[3];

    void toVertex(){
        for (int i=0; i<3; i++){
            Vertex v;
            v.setXYZ(X[i], Y[i], Z[i]);
            v.setColor(colors[i][0], colors[i][1], colors[i][2]);
            v.setNormal(normals[i][0], normals[i][1], normals[i][2]);
            v.setShading(shadings[i]);
            vertices[i] = v;
        }
    }

    void toBasic(){
        for (int i=0; i<3; i++){
            X[i] = vertices[i].X;
            Y[i] = vertices[i].Y;
            Z[i] = vertices[i].Z;
            colors[i][0] = vertices[i].color[0];
            colors[i][1] = vertices[i].color[1];
            colors[i][2] = vertices[i].color[2];
            normals[i][0] = vertices[i].normal[0];
            normals[i][1] = vertices[i].normal[1];
            normals[i][2] = vertices[i].normal[2];
            shadings[i] = vertices[i].shading;
        }
    }

    void shading2Basic(){
        for (int i=0; i<3; i++){
            shadings[i] = vertices[i].shading;
        }
    }

    void printTriangle(){
        cout<<std::setprecision(10);
        for (int i=0; i<3; i++){
            vertices[i].printVertex();
            cout<<"triangle print\n";
            cout<<X[i]<<", "<<Y[i]<<", "<<Z[i]<<"\n";
            cout<<"color: ";
            cout<<ceil(colors[i][0]*255)<<", ";
            cout<<ceil(colors[i][1] * 255)<<", "<<ceil(colors[i][2] * 255)<<"\n";
            cout<<"normals: ";
            cout<<normals[i][0]<<", "<<normals[i][1]<<", "<<normals[i][2]<<"\n";
            cout<<"Shading: "<<shadings[i]<<"\n\n";
        }
        
        cout<<"------------------------------\n\n";
    }

    void setCPosition(double c1, double c2, double c3){
        cPosition[0] = c1;
        cPosition[1] = c2;
        cPosition[2] = c3;
    }

    bool isGoingDown(){
        if (vertices[0].Y == vertices[1].Y && vertices[1].Y > vertices[2].Y){
            return true;
        }
        else if (vertices[1].Y == vertices[2].Y && vertices[1].Y > vertices[0].Y){
            return true;
        }
        else if (vertices[0].Y == vertices[2].Y && vertices[0].Y > vertices[1].Y){
            return true;
        }
        else{
            return false;
        }
    }

    bool isGoingUp(){
        if (vertices[0].Y == vertices[1].Y && vertices[1].Y < vertices[2].Y){
            return true;
        }
        else if (vertices[1].Y == vertices[2].Y && vertices[1].Y < vertices[0].Y){
            return true;
        }
        else if (vertices[0].Y == vertices[2].Y && vertices[0].Y < vertices[1].Y){
            return true;
        }
        else{
            return false;
        }
    }

    Vertex getFourthVertex(){

        Vertex v;

        double x1 = X[getLowVertexIndex()];
        double y1 = Y[getLowVertexIndex()];
        double z1 = Z[getLowVertexIndex()];
        double x2 = X[getHighVertexIndex()];
        double y2 = Y[getHighVertexIndex()];
        double z2 = Z[getHighVertexIndex()];
        double y = Y[getMiddleVertexIndex()];

        double x = getEndValue(x1, x2, y1, y2, y);
        double z = getEndValue(z1, z2, y1, y2, y);

        v.setXYZ(x,y,z);
        return interpolate(vertices[getLowVertexIndex()], vertices[getHighVertexIndex()] , v);
    
    }

    //Only for arbitrary
    int getLowVertexIndex(){
        return smallestIndex();
    }

    //Only for arbitrary
    int getHighVertexIndex(){
        return largestIndex();
    }

    int getMiddleVertexIndex(){
        int lowidx = getLowVertexIndex();
        int hihghidx = getHighVertexIndex();
        for (int i=0; i<3; i++){
            if (i != lowidx && i != hihghidx){
                return i;
            }
        }
        return -1;
    }

    //index of 2 triangles on the left
    int* get2LeftXIndex(){
        static int result[2];
        //lowest vertex
        result[0] = getExtremeIndex();

        int idx1 = getIndexOfRest()[0]; 
        int idx2 = getIndexOfRest()[1];
        if (vertices[idx1].X < vertices[idx2].X){
            //vertex on the left
            result[1] = idx1;
        }
        else{
            result[1] = idx2;
        }
        return result;
    }

    //index of 2 triangles on the right
    int* get2RightXIndex(){
        static int result[2];
        //lowest vertex
        result[0] = getExtremeIndex();

        int idx1 = getIndexOfRest()[0]; 
        int idx2 = getIndexOfRest()[1];
        if (vertices[idx1].X > vertices[idx2].X){
            //vertex on the right
            result[1] = idx1;
        }
        else{
            result[1] = idx2;
        }
        return result;
    }

  private:
    int getExtremeIndex(){
        int idx = 0;
        for (int i=0; i<3; i++){
            if(vertices[i].Y != vertices[(i+1)%3].Y 
                && vertices[i].Y != vertices[(i+2)%3].Y){
                idx = i;
            }
        }
        return idx;
    }

    int smallestIndex(){
        int idx = 0;
        for (int k=0; k<3; k++){
            if (vertices[k].Y < vertices[idx].Y){
                idx = k;
            }
        }
        return idx;
    }

    int largestIndex(){
        int idx = 0;
        for (int k=0; k<3; k++){
            if (vertices[k].Y > vertices[idx].Y){
                idx = k;
            }
        }
        return idx;
    }

    int* getIndexOfRest(){
        static int result[2];
        int k = getExtremeIndex();
        int j = 0;
        for (int i=0; i<3; i++){
            if (i != k){
                result[j] = i;
                j++;
            }
        }
        return result;
    }

    // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      double   *depth;
      int width, height;

      void InitializeScreen(){
        int npixels = width * height;
        for (int i = 0 ; i < npixels; i++){
            buffer[3 * i] = 0;
            buffer[3*i + 1] = 0;
            buffer[3*i + 2] = 0;
            depth[i] = -1.0;
        }
      }

  // would some methods for accessing and setting pixels be helpful?
};

Vertex getLeftEnd(Triangle t, int y){

    double x1 = t.X[t.get2LeftXIndex()[0]];
    double y1 = t.Y[t.get2LeftXIndex()[0]];
    double z1 = t.Z[t.get2LeftXIndex()[0]];
    double x2 = t.X[t.get2LeftXIndex()[1]];
    double y2 = t.Y[t.get2LeftXIndex()[1]];
    double z2 = t.Z[t.get2LeftXIndex()[1]];

    Vertex v;
    double x = getEndValue(x1, x2, y1, y2, y);
    double z = getEndValue(z1, z2, y1, y2, y);

    v.setXYZ(x, y, z);
    return interpolate(t.vertices[t.get2LeftXIndex()[0]], t.vertices[t.get2LeftXIndex()[1]], v);
}

Vertex getRightEnd(Triangle t, int y){

    double x1 = t.X[t.get2RightXIndex()[0]];
    double y1 = t.Y[t.get2RightXIndex()[0]];
    double z1 = t.Z[t.get2RightXIndex()[0]];
    double x2 = t.X[t.get2RightXIndex()[1]];
    double y2 = t.Y[t.get2RightXIndex()[1]];
    double z2 = t.Z[t.get2RightXIndex()[1]];

    Vertex v;
    double x = getEndValue(x1, x2, y1, y2, y);
    double z = getEndValue(z1, z2, y1, y2, y);

    v.setXYZ(x, y, z);

    return interpolate(t.vertices[t.get2RightXIndex()[0]], t.vertices[t.get2RightXIndex()[1]], v);
}

int RasterizeGoingDownTriangle(Triangle t, Screen s){

    double lowY = std::min({t.Y[0], t.Y[1], t.Y[2] });
    double highY = std::max({t.Y[0], t.Y[1], t.Y[2] });

    if (highY < 0){
        return -1;
    }
    if (lowY > s.height){
        lowY = s.height - 1;
    }

    for (int row = (int) ceil(lowY); row<= (int) floor(highY); row++){
        if (row < 0 || row >= s.height){
            continue;
        }
        int pixelRow = row * s.width * 3;
        Vertex leftEnd = getLeftEnd(t, row);
        Vertex rightEnd = getRightEnd(t, row);

        if (ceil(leftEnd.getX()) > floor(rightEnd.getX())){
            continue;
        }

        for (int col = ceil(leftEnd.getX()); col <=floor(rightEnd.getX()); col++){
            if (col >=0 && col <s.width){
                Vertex v;
                double z = getEndValue(leftEnd.getZ(), rightEnd.getZ(), leftEnd.getX(), rightEnd.getX(), col);
                v.setXYZ(col, row, z);

                Vertex vk = interpolate(leftEnd, rightEnd, v);
                
                if(vk.getZ() > s.depth[row * s.width + col]){

                    s.depth[row * s.width + col] = vk.getZ();

                    double tempR = std::min( 1.0, vk.color[0] * vk.shading );
                    double tempG = std::min( 1.0, vk.color[1] * vk.shading );
                    double tempB = std::min( 1.0, vk.color[2] * vk.shading );

                    s.buffer[pixelRow + 3*col] = ceil(255 * tempR );
                    s.buffer[pixelRow + 3*col + 1] = ceil(255 * tempG );
                    s.buffer[pixelRow + 3*col + 2] = ceil(255 * tempB );
                }   
            }
        }
    }
    return -1;
}

int RasterizeGoingUpTriangle(Triangle t, Screen s){

    double lowY = std::min({t.Y[0], t.Y[1], t.Y[2] });
    double highY = std::max({t.Y[0], t.Y[1], t.Y[2] });

    if (lowY > s.height){
        return -1;
    }
    if (highY < 0){
        highY = 0;
    }

    for (int row = ceil(lowY); row <= floor(highY); row++){
        if (row < 0 || row >= s.height){
            continue;
        }
        int pixelRow = row * s.width * 3;
        Vertex leftEnd = getLeftEnd(t, row);
        Vertex rightEnd = getRightEnd(t, row);

        if (ceil(leftEnd.getX()) > floor(rightEnd.getX())){
            continue;
        }

        for (int col = ceil(leftEnd.getX()); col <= floor(rightEnd.getX()); col++){
            if (col >=0 && col <s.width){
                Vertex v;
                double z = getEndValue(leftEnd.getZ(), rightEnd.getZ(), leftEnd.getX(), rightEnd.getX(), col);
                v.setXYZ(col, row, z);

                Vertex vk = interpolate(leftEnd, rightEnd, v);

                if(vk.getZ() > s.depth[row * s.width + col]){

                    s.depth[row * s.width + col] = vk.getZ();

                    double tempR = std::min( 1.0, vk.color[0] * vk.shading );
                    double tempG = std::min( 1.0, vk.color[1] * vk.shading );
                    double tempB = std::min( 1.0, vk.color[2] * vk.shading );

                    s.buffer[pixelRow + 3*col] = ceil(255 * tempR );
                    s.buffer[pixelRow + 3*col + 1] = ceil(255 * tempG );
                    s.buffer[pixelRow + 3*col + 2] = ceil(255 * tempB );
                }    
            }
        } 
    }
    return -1;
}

void RasterizeTriangle(Triangle t, Screen &s){
    
    if(t.isGoingDown()){
        RasterizeGoingDownTriangle(t, s);
    }
    else if(t.isGoingUp()){
        RasterizeGoingUpTriangle(t, s);
    }
    else{
        Triangle tDown;
        Triangle tUp;

        tDown.vertices[0] = t.vertices[t.getLowVertexIndex()];
        tDown.vertices[1] = t.vertices[t.getMiddleVertexIndex()];
        tDown.vertices[2] = t.getFourthVertex();
        tDown.shading2Basic();
        tDown.toBasic();


        tUp.vertices[0] = t.vertices[t.getHighVertexIndex()];
        tUp.vertices[1] = t.vertices[t.getMiddleVertexIndex()];
        tUp.vertices[2] = t.getFourthVertex();
        tUp.shading2Basic();
        tUp.toBasic();


        RasterizeGoingUpTriangle(tUp, s);
        RasterizeGoingDownTriangle(tDown, s);
    }
}

class Matrix
{
  public:
    double          A[4][4];

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;          // the closest distance from origin to each plan
    double          angle;              //
    double          position[3];        // position of camera
    double          focus[3];           // focual point between near and far 
    double          up[3];              // which direction is up going, perpendicular to the near vector

    double          U[3];
    double          V[3];
    double          W[3];
    double          T[3];

    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(int n, int m);

    void            config();
};

Matrix Camera::ViewTransform(){
    Matrix m;
    for (int i = 0 ; i < 4 ; i++){
        for (int j = 0 ; j < 4 ; j++)
        {
            m.A[i][j] = 0;
        }
    }
    m.A[0][0] = 1 / tan(angle / 2);
    m.A[1][1] = 1 / tan(angle / 2);
    m.A[2][2] = (far + near) / (far - near);
    m.A[2][3] = -1;
    m.A[3][2] = 2 * far * near / (far - near);

    return m;
}

Matrix Camera::CameraTransform(void) {
    Matrix m;

    double ut = 0;
    double vt = 0;
    double wt = 0;
    for (int i=0; i<3; i++){
        ut = ut + U[i] * T[i];
        vt = vt + V[i] * T[i];
        wt = wt + W[i] * T[i];
    }

    m.A[0][0] = U[0];
    m.A[1][0] = U[1];
    m.A[2][0] = U[2];
    m.A[3][0] = ut;

    m.A[0][1] = V[0];
    m.A[1][1] = V[1];
    m.A[2][1] = V[2];
    m.A[3][1] = vt;

    m.A[0][2] = W[0];
    m.A[1][2] = W[1];
    m.A[2][2] = W[2];
    m.A[3][2] = wt;

    m.A[0][3] = 0;
    m.A[1][3] = 0;
    m.A[2][3] = 0;
    m.A[3][3] = 1;

    return m;
}

Matrix Camera::DeviceTransform(int n, int m) {
    Matrix mat;

    for (int i = 0 ; i < 4 ; i++){
        for (int j = 0 ; j < 4 ; j++)
        {
            mat.A[i][j] = 0;
        }
    }
    mat.A[0][0] = n/2;
    mat.A[1][1] = m/2;
    mat.A[2][2] = 1;
    mat.A[3][0] = n/2;
    mat.A[3][1] = m/2;
    mat.A[3][3] = 1;

    return mat;
}

void Camera::config(){

    // u = up x (O-focus) => up x w
    // v = (O-focus) x u  => w x u
    // w = (O-focus)
    // O

    for (int i=0; i<3; i++){
        W[i] = position[i] - focus[i];
        T[i] = - position[i];
    }

    normalize(W, W, 3);

    U[0] = up[1] * W[2] - up[2] * W[1];
    U[1] = up[2] * W[0] - up[0] * W[2];
    U[2] = up[0] * W[1] - up[1] * W[0];
    normalize(U, U, 3);

    V[0] = W[1] * U[2] - W[2] * U[1];
    V[1] = W[2] * U[0] - W[0] * U[2];
    V[2] = W[0] * U[1] - W[1] * U[0];
    normalize(V, V, 3);

}

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

Camera GetCamera(int frame, int nframes)
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

std::vector<Triangle> GetTriangles(void)
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
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
#endif
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

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
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

int main (){

    vtkImageData *image = NewImage(1000, 1000);
    int npixels = 1000 * 1000;
    
    double *depth = new double [npixels];
    unsigned char *buffer = 
       (unsigned char *) image->GetScalarPointer(0,0,0);

    for (int i = 0 ; i < npixels; i++){
        buffer[3 * i] = 0;
        buffer[3*i + 1] = 0;
        buffer[3*i + 2] = 0;
        depth[i] = -1.0;
    }

    // Allocate Screen
    Screen screen;
    screen.buffer = buffer;
    screen.depth = depth;
    screen.width = 1000;
    screen.height = 1000;

    int size = 1000;
    
    for (int i=0; i<size; i++){

        screen.InitializeScreen();
        Camera c = GetCamera(i, 1000);
        c.config();

        Matrix ct = c.CameraTransform();
        Matrix vt = c.ViewTransform();
        Matrix dt = c.DeviceTransform(screen.width, screen.height);
        Matrix W2I = ct.ComposeMatrices(ct, vt);
        Matrix W2D = W2I.ComposeMatrices(W2I, dt);

        int j=0;
        std::vector<Triangle> triangles = GetTriangles();
        for (std::vector<Triangle>::iterator it = triangles.begin(); it!=triangles.end();++it){
            it->toVertex();
            it->setCPosition(c.position[0], c.position[1], c.position[2]);

            double viewDirection[3];
            for (int l = 0; l<3; l++){
                viewDirection[0] = c.position[0] - it->vertices[l].getX();
                viewDirection[1] = c.position[1] - it->vertices[l].getY();
                viewDirection[2] = c.position[2] - it->vertices[l].getZ();
                normalize(viewDirection, viewDirection, 3);
                it->vertices[l].CalculatePhongShading(lp, viewDirection);
            }

            // Transform Triangles
            for (int k=0; k<3; k++){
                double v[4];
                double transformed[4];

                v[0] = it->X[k];
                v[1] = it->Y[k];
                v[2] = it->Z[k];
                v[3] = 1;
                normalize(v, v, 4);

                W2D.TransformPoint(v, transformed);
                it->X[k] = transformed[0] / transformed[3];
                it->Y[k] = transformed[1] / transformed[3];
                it->Z[k] = transformed[2] / transformed[3];
            }
            it->shading2Basic();
            it->toVertex();

            RasterizeTriangle(*it, screen);
            j++;
        }

        char name[20];
        sprintf(name, "./images/frame%03d", i);
        WriteImage(image, name); 
    }

}
