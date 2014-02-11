////////////////////////////////////////////////////////////////////////
//
//   Harvard University
//   CS175 : Computer Graphics
//   Professor Steven Gortler
//
////////////////////////////////////////////////////////////////////////
/****************************************************** 
* Project:         CS 116B Homework #1
* File:            hw1.cpp 
* Purpose:         To experiment with various kinds of texture mapping techniques.
* Start date:      2/10/14 
* Programmers:      Zane Melcho, Jason Hungerford, Cesar Inarrea
* 
****************************************************** 
*/

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <stack>
#if __GNUG__
#   include <tr1/memory>
#endif

#include <GL/glew.h>
#ifdef __MAC__
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif


#include "cvec.h"
#include "matrix4.h"
#include "geometrymaker.h"
#include "ppm.h"
#include "glsupport.h"
#include "rigtform.h"

#define M_PI 3.1415926535897932384626433832795;
enum {ASPECT, FOV, ZAXIS};
enum {DIFFUSE, SOLID, SHINY, ANISOTROPY};


using namespace std;      // for string, vector, iostream, and other standard C++ stuff
using namespace tr1; // for shared_ptr

// G L O B A L S ///////////////////////////////////////////////////

// --------- IMPORTANT --------------------------------------------------------
// Before you start working on this assignment, set the following variable
// properly to indicate whether you want to use OpenGL 2.x with GLSL 1.0 or
// OpenGL 3.x+ with GLSL 1.3.
//
// Set g_Gl2Compatible = true to use GLSL 1.0 and g_Gl2Compatible = false to
// use GLSL 1.3. Make sure that your machine supports the version of GLSL you
// are using. In particular, on Mac OS X currently there is no way of using
// OpenGL 3.x with GLSL 1.3 when GLUT is used.
//
// If g_Gl2Compatible=true, shaders with -gl2 suffix will be loaded.
// If g_Gl2Compatible=false, shaders with -gl3 suffix will be loaded.
// To complete the assignment you only need to edit the shader files that get
// loaded
// ----------------------------------------------------------------------------
static const bool g_Gl2Compatible = false;

// Forward Declarations
struct RigidBody;
struct ShaderState;
static RigidBody* buildIcosahedron();
static const ShaderState& setupShader(int material);

static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static float g_frustNear = -0.1;    // near plane
static float g_frustFar = -50.0;    // far plane

static const float g_groundY = -5.0;      // y coordinate of the ground
static const float g_groundSize = 20.0;   // half the ground length

static int g_windowWidth = 512;
static int g_windowHeight = 512;
//static double g_aspect = 5.5;//g_windowWidth / g_windowHeight;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;
static float g_treeHeight = 15;
static float g_lastTreeX = 0;
static const int g_numOfObjects = 1; //Number of objects to be drawn
static bool isKeyboardActive = true;
static int mode = ASPECT;

struct ShaderState {
  GlProgram program;

  // Handles to uniform variables
  GLint h_uLight, h_uLight2;
  GLint h_uProjMatrix;
  GLint h_uModelViewMatrix;
  GLint h_uNormalMatrix;
  GLint h_uColor;

  // Handles to vertex attributes
  GLint h_aPosition;
  GLint h_aNormal;
  GLint h_aTangent;

  ShaderState(const char* vsfn, const char* fsfn) {
    readAndCompileShader(program, vsfn, fsfn);

    const GLuint h = program; // short hand

    // Retrieve handles to uniform variables
    h_uLight = safe_glGetUniformLocation(h, "uLight");
    h_uLight2 = safe_glGetUniformLocation(h, "uLight2");
    h_uProjMatrix = safe_glGetUniformLocation(h, "uProjMatrix");
    h_uModelViewMatrix = safe_glGetUniformLocation(h, "uModelViewMatrix");
    h_uNormalMatrix = safe_glGetUniformLocation(h, "uNormalMatrix");
    h_uColor = safe_glGetUniformLocation(h, "uColor");

    // Retrieve handles to vertex attributes
    h_aPosition = safe_glGetAttribLocation(h, "aPosition");
    h_aNormal = safe_glGetAttribLocation(h, "aNormal");
	 h_aTangent = safe_glGetAttribLocation(h, "aTangent");

    if (!g_Gl2Compatible)
      glBindFragDataLocation(h, 0, "fragColor");
    checkGlErrors();
  }

};

static const int g_numShaders = 4;
static const char * const g_shaderFiles[g_numShaders][2] = {
	{"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
	{"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"},
	{"./shaders/basic-gl3.vshader", "./shaders/shiny-gl3.fshader"},
	{"./shaders/basic-gl3.vshader", "./shaders/anisotropy-gl3.fshader"}
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
  {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"},
	{"./shaders/basic-gl2.vshader", "./shaders/shiny-gl2.fshader"},
	{"./shaders/basic-gl2.vshader", "./shaders/anisotropy-gl2.fshader"}
};
static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states

// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// A vertex with floating point position and normal
struct VertexPN {
  Cvec3f p, n, t;

  VertexPN() {}
  VertexPN(float x, float y, float z,
           float nx, float ny, float nz,
			  float tx, float ty, float tz)
    : p(x,y,z), n(nx, ny, nz), t(tx, ty, tz)
  {}

  // Define copy constructor and assignment operator from GenericVertex so we can
  // use make* functions from geometrymaker.h
  VertexPN(const GenericVertex& v) {
    *this = v;
  }

  VertexPN& operator = (const GenericVertex& v) {
    p = v.pos;
    n = v.normal;
	 t = v.tangent;
    return *this;
  }
};

struct Geometry {
  GlBufferObject vbo, ibo;
  int vboLen, iboLen;

  Geometry(VertexPN *vtx, unsigned short *idx, int vboLen, int iboLen) {
    this->vboLen = vboLen;
    this->iboLen = iboLen;

    // Now create the VBO and IBO
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPN) * vboLen, vtx, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
  }

  void draw(const ShaderState& curSS) {
    // Enable the attributes used by our shader
    safe_glEnableVertexAttribArray(curSS.h_aPosition);
    safe_glEnableVertexAttribArray(curSS.h_aNormal);
	 safe_glEnableVertexAttribArray(curSS.h_aTangent);

    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, p));
    safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, n));
	 safe_glVertexAttribPointer(curSS.h_aTangent, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, t));

    // bind ibo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

    // draw!
    glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);

    // Disable the attributes used by our shader
    safe_glDisableVertexAttribArray(curSS.h_aPosition);
    safe_glDisableVertexAttribArray(curSS.h_aNormal);
	 safe_glDisableVertexAttribArray(curSS.h_aTangent);
  }

	void draw(const ShaderState& curSS, Matrix4 MVM)
	{
		Matrix4 NMVM = normalMatrix(MVM);

		GLfloat glmatrix[16];
		MVM.writeToColumnMajorMatrix(glmatrix); // send MVM
		safe_glUniformMatrix4fv(curSS.h_uModelViewMatrix, glmatrix);

		NMVM.writeToColumnMajorMatrix(glmatrix); // send NMVM
		safe_glUniformMatrix4fv(curSS.h_uNormalMatrix, glmatrix);

		draw(curSS);
	}
};
/*-----------------------------------------------*/
struct RigidBody
{
	RigTForm rtf;
	Matrix4 scale;
	RigidBody **children;
	int numOfChildren;
	Cvec3 color;
	Geometry *geom;
	bool isVisible;
	bool isChildVisible;
	string name;
	int material;

	RigidBody()
	{
		rtf = RigTForm();
		scale = Matrix4();
		children = NULL;
		numOfChildren = 0;
		color = Cvec3(.5,.5,.5);
		geom = NULL;
		isVisible = true;
		isChildVisible = true;
		material = SOLID;
	}

	~RigidBody()
	{
		for (int i =0; i < numOfChildren; i++)
			delete children[i];
		delete []children;
		delete geom;
	}

	RigidBody(RigTForm rtf_, Matrix4 scale_, RigidBody **children_, Geometry *geom_, Cvec3 color_, int material_)
	{
		/* PURPOSE:		 
			RECEIVES:	 
							
			RETURNS:		
		*/

		rtf = rtf_;
		scale = scale_;
		children = children_;
		numOfChildren = 0;
		geom = geom_;
		color = color_;
		isVisible = true;
		material = material_;
	}

	void drawRigidBody(RigTForm invEyeRbt)
	{
		RigTForm respectFrame = invEyeRbt;
		draw(respectFrame, Matrix4());
	}

	void draw(RigTForm respectFrame_, Matrix4 respectScale_)
	{
		const ShaderState& curSS = setupShader(material);
		safe_glUniform3f(curSS.h_uColor, color[0], color[1], color[2]);
	
		// Draw Parent
		RigTForm respectFrame = respectFrame_ * rtf;
		Matrix4 respectScale = respectScale_ * scale;
		Matrix4 MVM = RigTForm::makeTRmatrix(respectFrame, respectScale);

		if (isVisible)
		{
			if (geom != NULL)
				geom->draw(curSS, MVM);
		}

		//Draw Children
		if (isChildVisible)
		{
			for (int i = 0; i < numOfChildren; i++)
			{
				children[i]->draw(respectFrame, respectScale);
			}
		}
		
	}

	void draw(Matrix4 respectFrame_)
	{
		const ShaderState& curSS = setupShader(material);
		safe_glUniform3f(curSS.h_uColor, color[0], color[1], color[2]);
			
		//Draw parent
		Matrix4 respectFrame = respectFrame_ * RigTForm::makeTRmatrix(rtf, scale);
		Matrix4 MVM = respectFrame;

		if (isVisible)
		{
			if (geom != NULL)
				geom->draw(curSS, MVM);
		}

		//Draw Children
		for (int i = 0; i < numOfChildren; i++)
		{
			children[i]->draw(respectFrame);
		}
	}
};
/*-----------------------------------------------*/
// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_sphere, g_triangle;

// --------- Scene

static Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2000, -3000.0, -5000.0);  // define two lights positions in world space
static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 0.0, 10.0)); // Default camera
static RigTForm g_eyeRbt = g_skyRbt; //Set the g_eyeRbt frame to be default as the sky frame
static RigidBody g_rigidBodies[g_numOfObjects]; // Array that holds each Rigid Body Object

///////////////// END OF G L O B A L S //////////////////////////////////////////////////
/*-----------------------------------------------*/
static Matrix4 lookAt(Cvec3f eyePosition, Cvec3f lookAtPosition, Cvec3f upVector)
{
	Cvec3f x, y, z, w;
	double m[16];

	//Different from the book but works correctly
	z = normalize(eyePosition - lookAtPosition);
	x = normalize(cross(upVector,z));
	y = cross(z,x);	

	int k = 0;

	for (int i = 0; i < 3; i++)
	{
		m[k] = x[i];
		k++;
		m[k] = y[i];
		k++;
		m[k] = z[i];
		k++;
		m[k] = eyePosition[i];
		k++;
	}

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;

	//return Matrix4();

	return Matrix4(m, true);
}
/*-----------------------------------------------*/
static float angleBetween(Cvec3 vectorOne, Cvec3 vectorTwo)
{
	float temp = dot(vectorOne, vectorTwo);
	float vOneNorm = norm(vectorOne);
	float vTwoNorm = norm(vectorTwo);
	temp = (temp / (vOneNorm * vTwoNorm));
	temp = acos(temp) * 180;
	temp /= M_PI;

	//cout << "angle = " << temp << "\n";
	//Matrix4::print(Matrix4::makeXRotation(temp));

	return temp;
}
/*-----------------------------------------------*/
static float lookAt(Cvec3 eyePosition, Cvec3 upPosition)
{
	return -(90 - angleBetween(eyePosition, upPosition));
}
/*-----------------------------------------------*/
static void lookAtOrigin()
{
	// Set angle to look at the origin
	Cvec3 eye = g_eyeRbt.getTranslation();
	Cvec3 up = Cvec3(0,1,0);
	g_eyeRbt.setRotation(Quat().makeXRotation(lookAt(eye,up)));
}
/*-----------------------------------------------*/
static void initCamera()
{
	Cvec3 eye = Cvec3(0.0, 0.0, 10.0);
	//Cvec3 at = Cvec3(0.0, 0.0, 0.0);
	//Cvec3 up = Cvec3(0.0,1.0,0.0);
	//g_skyRbt = lookAt(eye, at, up); // Default camera
	//g_skyRbt.setRotation(Quat().makeXRotation(lookAt(eye,up))); // TODO Change so lookat is done after conversion to Matrix
	g_skyRbt.setTranslation(eye);
	g_eyeRbt = g_skyRbt;

	// Initialize near and far
	//float z = -eye[2];
	//g_frustNear = z / 2.0;
	//g_frustFar = (3 * z) / 2.0;
}
/*-----------------------------------------------*/
static void initGround() {
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
  VertexPN vtx[4] = {
    VertexPN(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0, 0, 1, 0),
    VertexPN(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0, 0, 1, 0),
  };
  unsigned short idx[] = {0, 1, 2, 0, 2, 3};
  g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6));
}
/*-----------------------------------------------*/
static Geometry* initCubes() 
{
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  return new Geometry(&vtx[0], &idx[0], vbLen, ibLen);
}
/*-----------------------------------------------*/
static Geometry* initTriangles() 
{
  int ibLen = 3;
  int vbLen = 3;

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeTriangle(vtx.begin(), idx.begin());
  return new Geometry(&vtx[0], &idx[0], vbLen, ibLen);
}
/*-----------------------------------------------*/
static Geometry* initIcosahedron() 
{
  int ibLen = 60;
  int vbLen = 12;

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeIcosahedron(vtx.begin(), idx.begin());
  return new Geometry(&vtx[0], &idx[0], vbLen, ibLen);
}
/*-----------------------------------------------*/
static Geometry* initSpheres() 
{
	int slices = 20;
	int stacks = 20;
	float radius = 1;
	int ibLen, vbLen;
	getSphereVbIbLen(slices, stacks, vbLen, ibLen);

	// Temporary storage for cube geometry
	vector<VertexPN> vtx(vbLen);
	vector<unsigned short> idx(ibLen);

	makeSphere(radius, slices, stacks, vtx.begin(), idx.begin());
	return new Geometry(&vtx[0], &idx[0], vbLen, ibLen);
}
/*-----------------------------------------------*/
static Geometry* initCylinders()
{
	float radius = 1;
	float height = 1;
	int slices = 20;
	int ibLen, vbLen;
	getCylinderVbIbLen(slices, vbLen, ibLen);

	// Temporary storage for cube geometry
	vector<VertexPN> vtx(vbLen);
	vector<unsigned short> idx(ibLen);

	makeCylinder(slices, radius, height, vtx.begin(), idx.begin());
	return new Geometry(&vtx[0], &idx[0], vbLen, ibLen);
}
/*-----------------------------------------------*/
static void initDie()
{
	/* PURPOSE:		Creates a Icosahedron die with an image textured surface 
		REMARKS:		Icosahedron is edge length 2
	*/

	RigidBody *die;
	die = buildIcosahedron();
	g_rigidBodies[0] = *die;

	glutPostRedisplay();
}
/*-----------------------------------------------*/
static RigidBody* buildIcosahedron()
{
	/* PURPOSE:		Creates an Icosahedron object  
		RETURNS:    RigidBody* that points to the Icosahedron
	*/
	float width = 1;
	float height = 1;
	float thick = 1;

	RigTForm rigTemp = RigTForm(Cvec3(0, 0, 0));
	Matrix4 scaleTemp = Matrix4();
	
	// Make container
	RigidBody *icosahedron = new RigidBody(RigTForm(), Matrix4(), NULL, initCubes(), Cvec3(0.5, 0.5, 0.5), DIFFUSE);
	icosahedron->isVisible = false;
	icosahedron->name = "container";

	// Make Icosahedron
	rigTemp = RigTForm(Cvec3(0, 0, 0));
	scaleTemp = Matrix4::makeScale(Cvec3(width, height, thick));

	RigidBody *die = new RigidBody(rigTemp, scaleTemp, NULL, initIcosahedron(), Cvec3(1,0,0), SHINY);
	die->name = "icosahedron";

	//Setup Children
	icosahedron->numOfChildren = 1;
	icosahedron->children = new RigidBody*[icosahedron->numOfChildren];
	icosahedron->children[0] = die;

	return icosahedron;
}
/*-----------------------------------------------*/
// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
  GLfloat glmatrix[16];
  projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
  safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
}

// takes MVM and its normal matrix to the shaders
static void sendModelViewNormalMatrix(const ShaderState& curSS, const Matrix4& MVM, const Matrix4& NMVM) {
  GLfloat glmatrix[16];
  MVM.writeToColumnMajorMatrix(glmatrix); // send MVM
  safe_glUniformMatrix4fv(curSS.h_uModelViewMatrix, glmatrix);

  NMVM.writeToColumnMajorMatrix(glmatrix); // send NMVM
  safe_glUniformMatrix4fv(curSS.h_uNormalMatrix, glmatrix);
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
	  g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}
/*-----------------------------------------------*/
static void drawStuff() 
{
	/* PURPOSE:		Draws objects in relative 3d space  
	*/

	// short hand for current shader state
	const ShaderState& curSS = *g_shaderStates[g_activeShader];

	// build & send proj. matrix to vshader
	const Matrix4 projmat = makeProjectionMatrix();
	sendProjectionMatrix(curSS, projmat);

	// Use the g_eyeRbt as the eyeRbt;
	const RigTForm eyeRbt = g_eyeRbt;
	const RigTForm invEyeRbt = inv(eyeRbt);

	const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
	const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
	safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
	safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

	// draw ground
	// ===========
	//
	const RigTForm groundRbt = RigTForm();  // identity
	Matrix4 MVM = RigTForm::makeTRmatrix(invEyeRbt * groundRbt, Matrix4());
	Matrix4 NMVM = normalMatrix(MVM);
	sendModelViewNormalMatrix(curSS, MVM, NMVM);
	safe_glUniform3f(curSS.h_uColor, 0.1, 0.95, 0.1); // set color
	g_ground->draw(curSS);

	// Draw all Rigid body objects
	for (int i = 0; i < g_numOfObjects; i++)
	{
		g_rigidBodies[i].drawRigidBody(invEyeRbt);
	}
}
/*-----------------------------------------------*/
static const ShaderState& setupShader(int material)
{
	/* PURPOSE:		Setups Shader based on material 
		RECEIVES:   material - enum value of shader to be used 
		RETURNS:    curSS - ShaderState to be used to draw object 
	*/

	// Current Shader State
	glUseProgram(g_shaderStates[material]->program);
	const ShaderState& curSS = *g_shaderStates[material];

	// Build & send proj. matrix to vshader
	const Matrix4 projmat = makeProjectionMatrix();
	sendProjectionMatrix(curSS, projmat);

	// Use the g_eyeRbt as the eyeRbt;
	const RigTForm eyeRbt = g_eyeRbt;
	const RigTForm invEyeRbt = inv(eyeRbt);

	const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
	const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
	safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
	safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

	const RigTForm identityRbt = RigTForm();
	Matrix4 MVM = RigTForm::makeTRmatrix(invEyeRbt * identityRbt, Matrix4());
	Matrix4 NMVM = normalMatrix(MVM);

	sendModelViewNormalMatrix(curSS, MVM, NMVM);

	return curSS;
}
/*-----------------------------------------------*/
static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff();

  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  glutPostRedisplay();
}
/*-----------------------------------------------*/
static void motion(const int x, const int y) {
	const double dx = x - g_mouseClickX;
	const double dy = g_windowHeight - y - 1 - g_mouseClickY;

	RigTForm m;
	if (g_mouseLClickButton && !g_mouseRClickButton) { // left button down?
		m = g_eyeRbt * RigTForm(Quat().makeXRotation(dy)) * RigTForm(Quat().makeYRotation(-dx)) * inv(g_eyeRbt);
	}
  else if (g_mouseRClickButton && !g_mouseLClickButton) 
  { // right button down?
		m = g_eyeRbt * RigTForm(Cvec3(dx, dy, 0) * 0.01) * inv(g_eyeRbt); //Update based on Eye Frame
	//m = Matrix4::makeTranslation(Cvec3(dx, dy, 0) * 0.01);
  }
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
    //m = g_eyeRbt * Matrix4::makeTranslation(Cvec3(0, 0, -dy) * 0.01) * inv(g_eyeRbt); //Update based on Eye Frame
	  m = g_eyeRbt * RigTForm(Cvec3(0, 0,dy) * 0.01) * inv(g_eyeRbt); //Update based on Eye Frame
  }

  if (g_mouseClickDown) {
//	  g_objectRbt[0] *= m; // Simply right-multiply is WRONG
	  //g_rigidBodies[0].rtf = m * g_rigidBodies[0].rtf; //Update Icosahedron
	  g_eyeRbt = m * g_eyeRbt; //Update Camera
	  glutPostRedisplay(); // we always redraw if we changed the scene
  }

  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}
/*-----------------------------------------------*/
static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;
}
/*-----------------------------------------------*/
static void keyboard(const unsigned char key, const int x, const int y) 
{
	/* PURPOSE:		OpenGL Callback for Keyboard presses 
		RECEIVES:	unsigned char key - Key pressed
						int x - Mouse x_position when key pressed
						int y - Mouse y_position when key pressed
		REMARKS:		Handles robot modifications based on key presses and then requests a redisplay
	*/

	if (isKeyboardActive)
	{
		switch (key) 
		{
			case 27:
				exit(0);                                  // ESC
			case 'h':
				cout << " ============== H E L P ==============\n\n"
				<< "h\t\thelp menu\n"
				<< "s\t\tsave screenshot\n"
				<< "m\t\tToggle flat shading on/off.\n"
				<< "o\t\tCycle object to edit\n"
				<< "v\t\tCycle view\n"
				<< "drag left mouse to rotate\n" << endl;
				break;
			case 's':
				glFlush();
				writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
				break;
			case 'm':
				g_activeShader ^= 1;
				break;
	  }
	
		if (key == 'a')
		{
			mode = ASPECT;
		}
		else if (key == 'f')
		{
			mode = FOV;
		}
		else if (key == 'z')
		{
			mode = ZAXIS;
		}
		else if (key == 'b')
		{
		}
		else if (key == '-')
		{
			if (mode == ASPECT)
			{
				//g_aspect -= 0.1;
				//g_windowWidth = g_aspect * g_windowHeight;
			}
			else if (mode == FOV)
			{
				//g_frustFovY -= 0.1;
			}
			else if (mode == ZAXIS)
			{
				// Move camera along the positive z-axis
				Cvec3 cameraTrans = g_eyeRbt.getTranslation();
				g_eyeRbt.setTranslation(cameraTrans + Cvec3(0,0,1));
			}
		}
		else if (key == '=')
		{
			if (mode == ASPECT)
			{
				//g_aspect += 0.1;
				//g_windowWidth = g_aspect * g_windowHeight;
			}
			else if (mode == FOV)
			{
				//g_frustFovY += 0.1;
			}
			else if (mode == ZAXIS)
			{
				// Move camera along the negative z-axis
				Cvec3 cameraTrans = g_eyeRbt.getTranslation();
				g_eyeRbt.setTranslation(cameraTrans - Cvec3(0,0,1));
			}
		}
	}

	glutPostRedisplay();
}
/*-----------------------------------------------*/
static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("Assignment 2");                       // title the window

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
}

static void initGLState() {
  glClearColor(128./255., 200./255., 255./255., 0.);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
  if (!g_Gl2Compatible)
    glEnable(GL_FRAMEBUFFER_SRGB);
}
/*-----------------------------------------------*/
static void initShaders() {
  g_shaderStates.resize(g_numShaders);
  for (int i = 0; i < g_numShaders; ++i) {
    if (g_Gl2Compatible)
	 {
		g_shaderStates[i].reset(new ShaderState(g_shaderFilesGl2[i][0], g_shaderFilesGl2[i][1]));
	 }
	 else
      g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
  }
}
/*-----------------------------------------------*/
static void initGeometry() 
{
	//Initialize Object Matrix array
	initDie();
	initGround();
}
/*-----------------------------------------------*/
static void initLights()
{
	g_light1 = Cvec3(g_lastTreeX, g_treeHeight + 5.0, 0);
}
/*-----------------------------------------------*/
int main(int argc, char * argv[]) {
  try {
		initGlutState(argc,argv);

		glewInit(); // load the OpenGL extensions

		cout << (g_Gl2Compatible ? "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.3") << endl;
		if ((!g_Gl2Compatible) && !GLEW_VERSION_3_0)
		throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
		else if (g_Gl2Compatible && !GLEW_VERSION_2_0)
		throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.0");

		initGLState();
		initShaders();
		initCamera();
		initGeometry();
		initLights();

		glutMainLoop();
		return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    cin.ignore(); // pause to let users read whatever error came out
    return -1;
  }
}
