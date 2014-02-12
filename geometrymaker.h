#ifndef GEOMETRYMAKER_H
#define GEOMETRYMAKER_H

#include <cmath>

#include "cvec.h"

//--------------------------------------------------------------------------------
// Helpers for creating some special geometries such as plane, cubes, and spheres
//--------------------------------------------------------------------------------


// A generic vertex structure containing position, normal, and texture information
// Used by make* functions to pass vertex information to the caller
struct GenericVertex {
  Cvec3f pos;
  Cvec3f normal;
  Cvec2f tex;
  Cvec3f tangent, binormal;

  GenericVertex() {}
  GenericVertex(
    float x, float y, float z,
    float nx, float ny, float nz,
    float tu, float tv,
    float tx, float ty, float tz,
    float bx, float by, float bz)
    : pos(x,y,z), normal(nx,ny,nz), tex(tu, tv), tangent(tx, ty, tz), binormal(bx, by, bz){}

	GenericVertex(const GenericVertex& v)
	{
		*this = v;
	}

	GenericVertex& operator = (const GenericVertex& v)
	{
		pos = v.pos;
		normal = v.normal;
		tex = v.tex;
		tangent = v.tangent;
		binormal = v.binormal;
		return *this;
	}

	void setTex(float tu, float tv)
	{
		tex[0] = tu;
		tex[1] = tv;
	}
};

inline void getPlaneVbIbLen(int& vbLen, int& ibLen) {
  vbLen = 4;
  ibLen = 6;
}

template<typename VtxOutIter, typename IdxOutIter>
void makePlane(float size, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float h = size / 2.0;
  *vtxIter = GenericVertex(    -h, 0, -h, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex(-h, 0,  h, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( h, 0,  h, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( h, 0, -h, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, -1);
  *idxIter = 0;
  *(++idxIter) = 1;
  *(++idxIter) = 2;
  *(++idxIter) = 0;
  *(++idxIter) = 2;
  *(++idxIter) = 3;
}

inline void getCubeVbIbLen(int& vbLen, int& ibLen) {
  vbLen = 24;
  ibLen = 36;
}

template<typename VtxOutIter, typename IdxOutIter>
void makeCube(float size, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float h = size / 2.0;
#define DEFV(x, y, z, nx, ny, nz, tu, tv) { \
    *vtxIter = GenericVertex(x h, y h, z h, \
                             nx, ny, nz, tu, tv, \
                             tan[0], tan[1], tan[2], \
                             bin[0], bin[1], bin[2]); \
    ++vtxIter; \
}
  Cvec3f tan(0, 1, 0), bin(0, 0, 1);
  DEFV(+, -, -, 1, 0, 0, 0, 0); // facing +X
  DEFV(+, +, -, 1, 0, 0, 1, 0);
  DEFV(+, +, +, 1, 0, 0, 1, 1);
  DEFV(+, -, +, 1, 0, 0, 0, 1);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(0, 1, 0);
  DEFV(-, -, -, -1, 0, 0, 0, 0); // facing -X
  DEFV(-, -, +, -1, 0, 0, 1, 0);
  DEFV(-, +, +, -1, 0, 0, 1, 1);
  DEFV(-, +, -, -1, 0, 0, 0, 1);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(1, 0, 0);
  DEFV(-, +, -, 0, 1, 0, 0, 0); // facing +Y
  DEFV(-, +, +, 0, 1, 0, 1, 0);
  DEFV(+, +, +, 0, 1, 0, 1, 1);
  DEFV(+, +, -, 0, 1, 0, 0, 1);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 0, 1);
  DEFV(-, -, -, 0, -1, 0, 0, 0); // facing -Y
  DEFV(+, -, -, 0, -1, 0, 1, 0);
  DEFV(+, -, +, 0, -1, 0, 1, 1);
  DEFV(-, -, +, 0, -1, 0, 0, 1);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 1, 0);
  DEFV(-, -, +, 0, 0, 1, 0, 0); // facing +Z
  DEFV(+, -, +, 0, 0, 1, 1, 0);
  DEFV(+, +, +, 0, 0, 1, 1, 1);
  DEFV(-, +, +, 0, 0, 1, 0, 1);

  tan = Cvec3f(0, 1, 0);
  bin = Cvec3f(1, 0, 0);
  DEFV(-, -, -, 0, 0, -1, 0, 0); // facing -Z
  DEFV(-, +, -, 0, 0, -1, 1, 0);
  DEFV(+, +, -, 0, 0, -1, 1, 1);
  DEFV(+, -, -, 0, 0, -1, 0, 1);
#undef DEFV

  for (int v = 0; v < 24; v +=4) {
    *idxIter = v;
    *++idxIter = v + 1;
    *++idxIter = v + 2;
    *++idxIter = v;
    *++idxIter = v + 2;
    *++idxIter = v + 3;
    ++idxIter;
  }
}

inline void getSphereVbIbLen(int slices, int stacks, int& vbLen, int& ibLen) {
  assert(slices > 1);
  assert(stacks >= 2);
  vbLen = (slices + 1) * (stacks + 1);
  ibLen = slices * stacks * 6;
}

inline void getCylinderVbIbLen(int slices, int& vbLen, int& ibLen) {
  assert(slices > 1);
  vbLen = slices * 2;
  ibLen = slices * 6; 
}

template<typename VtxOutIter, typename IdxOutIter>
void makeSphere(float radius, int slices, int stacks, VtxOutIter vtxIter, IdxOutIter idxIter) {
  using namespace std;
  assert(slices > 1);
  assert(stacks >= 2);

  const double radPerSlice = 2 * CS175_PI / slices;
  const double radPerStack = CS175_PI / stacks;

  vector<double> longSin(slices+1), longCos(slices+1);
  vector<double> latSin(stacks+1), latCos(stacks+1);
  for (int i = 0; i < slices + 1; ++i) {
    longSin[i] = sin(radPerSlice * i);
    longCos[i] = cos(radPerSlice * i);
  }
  for (int i = 0; i < stacks + 1; ++i) {
    latSin[i] = sin(radPerStack * i);
    latCos[i] = cos(radPerStack * i);
  }

  for (int i = 0; i < slices + 1; ++i) {
    for (int j = 0; j < stacks + 1; ++j) {
      float x = longCos[i] * latSin[j];
      float y = longSin[i] * latSin[j];
      float z = latCos[j];

      Cvec3f n(x, y, z);
      Cvec3f t(-longSin[i], longCos[i], 0);
      Cvec3f b = cross(n, t);

      *vtxIter = GenericVertex(
        x * radius, y * radius, z * radius,
        x, y, z,
        1.0/slices*i, 1.0/stacks*j,
        t[0], t[1], t[2],
        b[0], b[1], b[2]);
      ++vtxIter;

      if (i < slices && j < stacks ) {
        *idxIter = (stacks+1) * i + j;
        *++idxIter = (stacks+1) * i + j + 1;
        *++idxIter = (stacks+1) * (i + 1) + j + 1;

        *++idxIter = (stacks+1) * i + j;
        *++idxIter = (stacks+1) * (i + 1) + j + 1;
        *++idxIter = (stacks+1) * (i + 1) + j;
        ++idxIter;
      }
    }
  }
}

template<typename VtxOutIter, typename IdxOutIter>
void makeCylinder(int slices, float radius, float height, VtxOutIter vtxIter, IdxOutIter idxIter) 
{
	/* PURPOSE:		Fills a vertex and index buffer to make a cylinder 
		RECEIVES:   slices - number of points around the circle
						radius - radius of cylinder
						height - Height of cylinder
						vtxIter - VertexBuffer Iterator
						idxIter - indexBuffer Iterator
		REMARKS:    Creates a hollow cylinder 
	*/
  using namespace std;
  assert(slices > 1);
  int stacks = 2;  

  float angle = 2 * CS175_PI / slices;
  float halfHeight = height / 2;

  for (int i = 0; i < slices; ++i) 
  {
    for (int j = 0; j < stacks; ++j) 
	 { 
		float x = cos(i * angle);
		float z = sin(i * angle);
      
		int flip = -1;
		if (j == 0)
			flip *= -1;

		float y = halfHeight * flip;

      Cvec3f n(x, y, z);
      Cvec3f t(-x, y, 0);
      Cvec3f b = cross(n, t);

      *vtxIter = GenericVertex(
        x * radius, y, z * radius,
        x, y, z,
        1.0/slices*i, 1.0/2*j,
        t[0], t[1], t[2],
        b[0], b[1], b[2]);
      ++vtxIter;
    }
  }

  // Load the index buffer
	int k = 0;
	for (int i = 0; i < slices; i++)
	{
		int a = k % (slices * 2);
		int b = (k + 1) % (slices * 2);
		int c = (k + 2) % (slices * 2);
		int d = (k + 3) % (slices * 2);

		*idxIter = a;
		*++idxIter = c;
		*++idxIter = b;

		*++idxIter = b;
		*++idxIter = c;
		*++idxIter = d;
		++idxIter;
		k += 2;
	}
}

template<typename VtxOutIter, typename IdxOutIter>
void makeTriangle(VtxOutIter vtxIter, IdxOutIter idxIter) 
{
	using namespace std;

	float startX = sqrt(0.75) / 2.0;

	Cvec3f *threePoints = new Cvec3f[3];
	threePoints[0] = Cvec3f(-startX, 0, 0); 
	threePoints[1] = Cvec3f(startX, -0.5, 0);
	threePoints[2] = Cvec3f(startX, 0.5, 0);
	

	float textures[][2] = {{1,0},{0,0},{.5,1}};

	for (int i = 0; i < 3; i++)
	{
		float x = threePoints[i][0];
		float y = threePoints[i][1];
		float z = threePoints[i][2];

		Cvec3f n(x, y, z);
		Cvec3f t(-x, y, 0);
		Cvec3f b = cross(n, t);

		*vtxIter = GenericVertex(
			x, y, z,
			n[0], n[1], n[2],
			textures[i][0], textures[i][1],
			t[0], t[1], t[2],
			b[0], b[1], b[2]);
		++vtxIter;
	}
	
	idxIter[0] = 0;
	idxIter[1] = 1;
	idxIter[2] = 2;
}

template<typename VtxOutIter, typename IdxOutIter>
void makeIcosahedron(VtxOutIter vtxIter, IdxOutIter idxIter) 
{
	using namespace std;

	float phi = (1 + sqrt(5.0)) / 2.0;

	Cvec3f *points = new Cvec3f[12];

	points[0] = Cvec3f(-1, phi, 0);		//Top-left
	points[1] = Cvec3f(1, phi, 0);		//Top-right
	points[2] = Cvec3f(-1, -phi, 0);		//Bottom-left
	points[3] = Cvec3f(1, -phi, 0);		//Bottom-right
	
	points[4] = Cvec3f(0, 1, phi);		//Top-front
	points[5] = Cvec3f(0, 1, -phi);		//Top-back
	points[6] = Cvec3f(0, -1, phi);		//Bottom-front
	points[7] = Cvec3f(0, -1, -phi);		//Bottom-back

	points[8] = Cvec3f(-phi, 0, 1);		//Front-left
	points[9] = Cvec3f(phi, 0, 1);		//Front-right 
	points[10] = Cvec3f(-phi, 0, -1);	//Back-left
	points[11] = Cvec3f(phi, 0, -1);		//Back-right
	
	//Setup tu, tv
	float textures[][2] = {{0.5,1},{1,0},{0.5,1},{1,0},{0,0},{0,0},{0,0},{0,0},{1,0},{0.5,1},{0,0},{0,0}};
	//float textures[][2] = {{0,0},{1,1},{0,0},{0,0},{0,1},{0,0},{0,0},{0,0},{0,0},{1,0},{0,0},{0,0}};

	//Load all points into vtxIter
	for (int i = 0; i < 12; i++)
	{
		float x = points[i][0];
		float y = points[i][1];
		float z = points[i][2];

		Cvec3f n(x, y, z);
		Cvec3f t(-x, y, 0);
		Cvec3f b = cross(n, t);

		*vtxIter = GenericVertex(
			x, y, z,
			n[0], n[1], n[2],
			textures[i][0], textures[i][1],
			t[0], t[1], t[2],
			b[0], b[1], b[2]);
		++vtxIter;
	}

	int temp[] = {	0,10,8, //Top
						0,8,4,
						0,4,1,
						0,1,5,
						0,5,10,

						8,10,2, //Middle
						8,2,6,
						8,6,4,
						9,4,6,
						9,1,4,
						9,11,1,
						5,1,11,
						5,11,7,
						5,7,10,
						10,7,2,

						3,6,2, //Bottom
						3,9,6,
						3,11,9,
						3,7,11,
						3,2,7
					 };

	// Isolating one triangle
	int temp2[] = {0,6,9, //Top
						1,1,1,
						1,1,1,
						1,1,1,
						1,1,1,

						1,1,1, //Middle
						1,1,1,
						1,1,1,
						1,1,1,
						1,1,1,
						1,1,1,
						1,1,1,
						1,1,1,
						1,1,1,
						1,1,1,

						1,1,1, //Bottom
						1,1,1,
						1,1,1,
						1,1,1,
						1,1,1
					 };

	for (int i = 0; i < 60; i++)
	{
		idxIter[i] = temp[i];
	}
}

#endif
