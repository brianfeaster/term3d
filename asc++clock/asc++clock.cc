/* :NOTES:
 *  ctime
 */
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <unistd.h>
#include <ctime>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/ioctl.h>
#include <stdio.h>

#define PI   3.14159265359

void SWAP (int &a, int &b) { b^=a^=b,a^=b; }

/* Keep track of
 */
template <class T>
class V : public std::vector<T> {
  public:
  V<T> () = default;
  V<T> (const std::vector<T> &that) : std::vector<T>(that) {}
  V<T> (const V<T> &that) : std::vector<T>(that) {}

  friend std::ostream& operator<<(std::ostream& o, const V& self) {
    for (int i=0; (i < self.size()); i+=4) {
      o << "{" << self[0+i]  << " " << self[1+i]  << " " << self[2+i]  << " " << self[3+i]  << "}" << std::endl;
    }
    return o;
  }

};



template <class T>
class M : public std::vector<T> {

  public:
  M<T> () = default;
  M<T> (const std::vector<T> &that) : std::vector<T>(that) {}
  M<T> (const M<T> &that) : std::vector<T>(that) {}

  friend std::ostream& operator<<(std::ostream& o, const M& self) {
    o << "[" << self[0]  << " " << self[1]  << " " << self[2]  << " " << self[3]  << "]" << std::endl;
    o << "[" << self[4]  << " " << self[5]  << " " << self[6]  << " " << self[7]  << "]" << std::endl;
    o << "[" << self[8]  << " " << self[9]  << " " << self[10] << " " << self[11] << "]" << std::endl;
    o << "[" << self[12] << " " << self[13] << " " << self[14] << " " << self[15] << "]" << std::endl << std::endl;
    return o;
  }

  M<T> &operator*(const T &c) {
    //std::transform(std::begin(this), std::end(this), std::begin(this), [c](T e){e*c;});
    for (auto i = this->begin(); i != this->end(); ++i) {
      *i *= c;
    }
    return *this;
  }

  V<T> operator*(const V<T> &v) {
    T a, b, c, d;
    V<T> u;
    u.resize(v.size());
    for (int i=0; (i < v.size()); i+=4) {
      a=v[i+0]; b=v[i+1]; c=v[i+2]; d=v[i+3];
      u[i+0] = (*this)[0]*a  + (*this)[1]*b  + (*this)[2]*c  + (*this)[3]*d;
      u[i+1] = (*this)[4]*a  + (*this)[5]*b  + (*this)[6]*c  + (*this)[7]*d;
      u[i+2] = (*this)[8]*a  + (*this)[9]*b  + (*this)[10]*c + (*this)[11]*d;
      u[i+3] = (*this)[12]*a + (*this)[13]*b + (*this)[14]*c + (*this)[15]*d;
    }
    return u;
  }


  V<T> operator*(const V<T> &v) const {
    T a, b, c, d;
    V<T> u;
    u.resize(v.size());
    for (int i=0; (i < v.size()); i+=4) {
      a=v[i+0]; b=v[i+1]; c=v[i+2]; d=v[i+3];
      u[i+0] = (*this)[0]*a  + (*this)[1]*b  + (*this)[2]*c  + (*this)[3]*d;
      u[i+1] = (*this)[4]*a  + (*this)[5]*b  + (*this)[6]*c  + (*this)[7]*d;
      u[i+2] = (*this)[8]*a  + (*this)[9]*b  + (*this)[10]*c + (*this)[11]*d;
      u[i+3] = (*this)[12]*a + (*this)[13]*b + (*this)[14]*c + (*this)[15]*d;
    }
    return u;
  }


  M<T> &operator*(const M<T> &m) {
    M<T> &self=*this;
    T a=self[0], b=self[1], c=self[2], d=self[3];
    self[0] = a*m[0] + b*m[4] + c*m[8] + d*m[12];
     self[1] = a*m[1] + b*m[5] + c*m[9] + d*m[13];
      self[2] = a*m[2] + b*m[6] + c*m[10] + d*m[14];
       self[3] = a*m[3] + b*m[7] + c*m[11] + d*m[15];
    a=self[4]; b=self[5]; c=self[6]; d=self[7];
    self[4] = a*m[0] + b*m[4] + c*m[8] + d*m[12];
     self[5] = a*m[1] + b*m[5] + c*m[9] + d*m[13];
      self[6] = a*m[2] + b*m[6] + c*m[10] + d*m[14];
       self[7] = a*m[3] + b*m[7] + c*m[11] + d*m[15];
    a=self[8]; b=self[9]; c=self[10]; d=self[11];
    self[8] = a*m[0] + b*m[4] + c*m[8] + d*m[12];
     self[9] = a*m[1] + b*m[5] + c*m[9] + d*m[13];
      self[10] = a*m[2] + b*m[6] + c*m[10] + d*m[14];
       self[11] = a*m[3] + b*m[7] + c*m[11] + d*m[15];
    a=self[12]; b=self[13]; c=self[14]; d=self[15];
    self[12] = a*m[0] + b*m[4] + c*m[8] + d*m[12];
     self[13] = a*m[1] + b*m[5] + c*m[9] + d*m[13];
      self[14] = a*m[2] + b*m[6] + c*m[10] + d*m[14];
       self[15] = a*m[3] + b*m[7] + c*m[11] + d*m[15];
    return self;
  }
    

  static M<T> identity (void);
  static M<T> rotx (double r);
  static M<T> roty (double r);
  static M<T> rotz (double r);
  static M<T> translate (double x, double y, double z);
  static M<T> scale (double x, double y, double z);
};

// Factories
template<class T> M<T> M<T>::identity(void) { return M<T>({1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0}); }
template<class T> M<T> M<T>::rotx(double r) { return M<T>({1.0,0.0,0.0,0.0, 0.0,cos(r),sin(r),0.0, 0.0,-sin(r),cos(r),0.0, 0.0,0.0,0.0,1.0}); }
template<class T> M<T> M<T>::roty(double r) { return M<T>({cos(r),0.0,-sin(r),0.0, 0.0,1.0,0.0,0.0, sin(r),0.0,cos(r),0.0, 0.0,0.0,0.0,1.0}); }
template<class T> M<T> M<T>::rotz(double r) { return M<T>({cos(r),-sin(r),0.0,0.0, sin(r),cos(r),0.0,0.0,  0.0,0.0,1.0,0.0,  0.0,0.0,0.0,1.0}); }
template<class T> M<T> M<T>::translate(double x, double y, double z) { return M<T>({1,0,0,x,  0,1,0,y,  0,0,1,z,  0,0,0,1}); }
template<class T> M<T> M<T>::scale(double x, double y, double z) { return M<T>({x,0,0,0,  0,y,0,0,  0,0,z,0,  0,0,0,1}); }

using mat=class M<double>;
using vec=class V<double>;



class Face {

  public:
  int  type;
  int  color;
  char ch;
  vec vertices;

  Face () : type(0), color(7), ch('.'), vertices() {}
  Face (int _type, int _color, char _ch, const mat &transform, const vec &face)
    : type(_type), color(_color), ch(_ch), vertices({-1,1,0,1, 1,1,0,1,  -1,-1,0,1}) {
    vertices.resize(face.size()+12);
    for (int i=0; (i < face.size()); ++i) {
      vertices[i+12] = face[i];
    }
    vertices = transform * vertices; // This requies operator* be const
  }

};



/* Display class
*/
class Display {
  struct winsize win;
  int rows2;
  int cols2;
  std::vector<Face> faces;
  std::vector<int> ranges;
  vec normals;
  vec vertices;
  int l=1;
public:
  Display () {
    probeUserTerminal();
  }
  void probeUserTerminal () { ioctl(1, TIOCGWINSZ, &win); rows2=win.ws_row/2; cols2=win.ws_col/2; }
  void renderLine(int y1, int x2, int y3, int x4, char ch);
  void renderPoly(int from, int to, vec u, char ch, int color) {
    double a,b,c,x,y,z,m,n,o;
    x=a=u[0+from]; y=b=u[1+from]; z=c=u[2+from];
    std::cout << "\e[0;3" << color << "m";
    for (int i=4+from; i<to; i+=4) {
      o=u[i+2]; n=u[i+1]; m=u[i+0];
      renderLine(y/(z)+0.5, x/(z)+0.5, n/(o)+0.5, m/(o)+0.5, ch);
      x=m; y=n; z=o;
    }
    renderLine(y/(z)+0.5, x/(z)+0.5, b/(c)+0.5, a/(c)+0.5, ch);
  }
  void renderPoints(int from, int to, vec u, char ch, int color) {
    double a,b,c;
    std::cout << "\e[4" << color << "m";
    for (int i=from; i<to; i+=4) {
      a=u[0+i]; b=u[1+i]; c=u[2+i];
      renderLine(b/(c)+0.5, a/(c)+0.5, b/(c)+0.5, a/(c)+0.5, ch);
    }
    std::cout << "\e[0m";
  }
  void renderLines(int from, int to, vec u, char ch, int color) {
    double a,b,c,x,y,z;
    std::cout << "\e[3" << color << "m";
    for (int i=from; i<to; i+=8) {
      a=u[0+from]; b=u[1+from]; c=u[2+from];
      x=u[4+from]; y=u[5+from]; z=u[6+from];
      renderLine(b/(c)+0.5, a/(c)+0.5, y/(z)+0.5, x/(z)+0.5, ch);
    }
  }
  void clearFaces () {
    faces.resize(0);
    ranges.resize(0);
    normals.resize(0);
    vertices.resize(0);
  }
  void addFace (const Face &face) {
    int i;
    faces.push_back(face);
    for (i=0; i<12; ++i) normals.push_back(face.vertices[i]);
    ranges.push_back(vertices.size()); // Keep track of first index
    for (i=12; i<face.vertices.size(); ++i) vertices.push_back(face.vertices[i]);
    ranges.push_back(vertices.size()); // Keep track of last index, exclusive
    //std::cout << normals;
    //std::cout << vertices;
  }
  void addRandom () {
    faces.push_back(Face());
    for (int i=0; i<12; ++i) normals.push_back(0); // ignored
    ranges.push_back(vertices.size()); // Keep track of first index
    for (int i=0; i<50; ++i) {
      vertices.push_back((random()%50-25)/20.0);
      vertices.push_back((random()%50-25)/20.0);
      vertices.push_back((random()%50-25)/20.0);
      vertices.push_back(1);
    }
    ranges.push_back(vertices.size()); // Keep track of last index, exclusive
  }
  void addMarquee () {
    int ch, colors[]={1,3,2,4,6,5}, c=0;
    double x=0, y=-1;
    while ((ch=fgetc(stdin)) != EOF) {
      if (ch == '\n') {
        x=0.0;
        y=y+.1;
        ++c;
      } else {
        addFace(Face(0, colors[c/4%6], ch, M<double>::identity(), {{x,y,0,1}}));
        x=x+.1;
      }
    }
  }

  bool facing (int idx, const vec &norms) {
    // Check Z component of the normal of the above two vectors
    //std::cout << "\e[H" << norms;
    idx *= 12;
    double d;

    d = norms[idx+3];
    double x0=norms[idx+0]/d;
    double y0=norms[idx+1]/d;
    double z0=norms[idx+2]/d;

    d = norms[idx+7];
    double x1=norms[idx+4]/d;
    double y1=norms[idx+5]/d;
    double z1=norms[idx+6]/d;

    d = norms[idx+11];
    double x2=norms[idx+8]/d;
    double y2=norms[idx+9]/d;
    double z2=norms[idx+10]/d;

    double va0 = x1/z1 - x0/z0;
    double va1 = y1/z1 - y0/z0;
    double vb0 = x2/z2 - x0/z0;
    double vb1 = y2/z2 - y0/z0;
    return true || (va0 * vb1 - va1 * vb0) <= 0.0;
  }
  void renderAll (void) {
    /* Make a face {given: xform} {add: normal culling vertices} {given: 1,1, 1,-1, -1,-1, -1,1}
       Apply the xform to the vertices.  Each vertex set saved separately.  culling vertices all xform first and
       used to determine if the given vertices should be xformed and rendered.
    */
      double p=30.0;
      probeUserTerminal();
      //                  perspective transformation.
      mat m = mat({1,0,0,0,  0,1,0,0,  0,0,1/p,1, 0,0,0,1})  * M<double>::scale(10,10,10) * M<double>::roty(-l/30.0) * M<double>::rotx(l/40.0);
      vec u = m * vertices;
      vec n = m * normals;
      std::cout << "\e[2J";
      for (int r=0; (r < ranges.size()); r+=2) {
        if (faces[r/2].type == 0) {
          renderPoints(ranges[r], ranges[r+1], u, faces[r/2].ch, faces[r/2].color);
        } else {
          if (facing(r/2, n)) { 
            if (faces[r/2].type==1) renderLines(ranges[r], ranges[r+1], u, faces[r/2].ch, faces[r/2].color);
            if (faces[r/2].type==2) renderPoly(ranges[r], ranges[r+1], u, faces[r/2].ch, faces[r/2].color);
          }
        }
      }
      ++l;
      usleep(1000*20);
  }
};

mat rotyPI2   = M<double>::roty(PI/2.0);
mat rotyPI    = M<double>::roty(PI);
mat roty3PI2  = M<double>::rotx(3.0*PI/2.0);

mat rotxPI2   = M<double>::rotx(PI/2.0);
mat rotx3PI2  = M<double>::rotx(3.0*PI/2.0);
mat transzn1 = M<double>::translate(0.0, 0.0, -1.0);

int main (void) {
  vec unitFrame {{1,1,0,1, 1,-1,0,1, -1,-1,0,1, -1,1,0,1}};
  Display display;
    display.clearFaces();
    display.addRandom();
    //display.addMarquee();
    display.addFace(Face(2, 2, 'a',                             transzn1, unitFrame));
    display.addFace(Face(2, 3, 'b',  rotyPI2                  * transzn1, unitFrame));
    display.addFace(Face(2, 4, 'c',  rotyPI                   * transzn1, unitFrame));
    display.addFace(Face(2, 5, 'd',  roty3PI2                 * transzn1, unitFrame));
    display.addFace(Face(2, 7, '1',  rotxPI2  * transzn1, unitFrame));
    display.addFace(Face(2, 7, '2',  rotx3PI2 * transzn1, unitFrame));
  for (;;) {
    //struct timeval tv;
    //gettimeofday(&tv, NULL);
    //double sec  = 2.0 * PI * (tv.tv_sec%60) / 60.0;
    //double min  = 2.0 * PI * ((tv.tv_sec/60)%60) / 60.0;
    //double hour = 2.0 * PI * ((tv.tv_sec/3600)%12) / 60.0;

    //display.addFace(Face(1, 1, '@',  M<double>::rotz(sec)     * M<double>::translate(0, 0 , -1.2), {{0,-.7,0,1, 0,-1,  0,1}}));
    //display.addFace(Face(1, 2, '@',  M<double>::rotz(min)     * M<double>::translate(0, 0 , -1.2), {{0,0,  0,1, 0,-0.5,0,1}}));
    //display.addFace(Face(1, 4, '@',  M<double>::rotz(hour)    * M<double>::translate(0, 0 , -1.2), {{0,0,  0,1, 0,-1,  0,1}}));
    //display.addFace(Face(2, 2, 'a',  M<double>::roty(0)       * M<double>::translate(0, 0 , -1), {{1,1,0,1, 1,-1,0,1, -1,-1,0,1, -1,1,0,1}}));
    //display.addFace(Face(2, 7, 'b',  M<double>::roty(PI/2)    * M<double>::translate(0, 0 , -1), {{1,1,0,1, 1,-1,0,1, -1,-1,0,1, -1,1,0,1}}));
    //display.addFace(Face(2, 7, 'c',  M<double>::roty(PI)      * M<double>::translate(0, 0 , -1), {{1,1,0,1, 1,-1,0,1, -1,-1,0,1, -1,1,0,1}}));
    //display.addFace(Face(2, 7, 'd',  M<double>::roty(PI+PI/2) * M<double>::translate(0, 0 , -1), {{1,1,0,1, 1,-1,0,1, -1,-1,0,1, -1,1,0,1}}));
    //display.addFace(Face(2, 7, '1',  M<double>::rotx(PI/2)    * M<double>::translate(0, 0 , -1), {{1,1,0,1, 1,-1,0,1, -1,-1,0,1, -1,1,0,1}}));
    //display.addFace(Face(2, 7, '2',  M<double>::rotx(PI+PI/2) * M<double>::translate(0, 0 , -1), {{1,1,0,1, 1,-1,0,1, -1,-1,0,1, -1,1,0,1}}));
    display.renderAll();
  }
  return 0;
}


void Display::renderLine (int y1, int x2, int y3, int x4, char chr) {
  int y, x, n, m, dy, dx, yy, yyxx, e, i, twx, twy, twsx, twsy;
  const char *w, *ws;

  // Make sure line is drawn from left -x to right +x so r eorder points so x < m
  if (x2 <= x4) {
    y=y1; x=x2; n=y3; m=x4;
  } else {
    y=y3; x=x4; n=y1; m=x2;
  }

  dx = m - x;
  dy = (y1 < y3) ? y3-y1 : y1-y3;

  if ( dy < dx ) { // Small slope...
    w = ""; // so walk right...
    ws = (y < n) ? "\v" : "\eM"; // and step down or up.

    twx=1;  twy=0;
    twsx=1;  twsy=(y < n) ? 1 : -1; // and step down or up.
  } else { // Large slope...
    if ( y < n ) { // so walk down and step right
       w = "\b\v";  twx=0; twy=1;
       ws = "\v";   twsx=1; twsy=1;
    } else { // or walk up and step right
       w = "\b\eM";  twx=0; twy=-1;
       ws = "\eM";   twsx=1; twsy=-1;
    }
    SWAP(dy, dx);
  }

  int ty=y+rows2;
  int tx=x+cols2;
  std::cout << "\e[" << ty << ";" << tx << "H" << chr;

  yy = 2*dy;
  yyxx = yy - 2*dx;
  e = yy - dx;
  i = dx;
  while (i--) {
    if (0 <= e) {
      std::cout << ws << chr;
      e = e + yyxx;
    } else {
      std::cout << w << chr;
      e = e + yy;
    }
  }

  std::cout << std::flush;
}
