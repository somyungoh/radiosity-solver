#define _CRT_SECURE_NO_WARNINGS

//**************************************************************************
//
//		Somyung(David) Oh
//
//		Project3 - Radiosity
//
//		Computer Graphics, 2017 Fall
//		Texas A&M University
//
//		<Projet Tasks>
//		1. Compute patch to element formfactor using Hemicube
//		2. Solve radiosity equation by progressive refinement approach
//		3. Compute ambient term
//		4. Interpolate vertex colors by neighboring elements
//
//**************************************************************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <time.h>
#include "math.h"
#include "GL\glui.h"
#include "glm\glm.hpp"
#include "GL\glut.h"

using namespace std;
using namespace glm;

#define HEMICUBE_HEIGHT 1.0
#define HEMICUBE_SUBDIV 512
#define PI 3.141592

typedef vec3 Color;
typedef vec3 Vertex;
typedef vec3 Vector;
enum { FRONT, LEFT, RIGHT, TOP, BOTTOM };


//		Global Variables		//

//	GLUI Variables
GLUI			 *glui;
GLUI_Panel		 *panel_control;
GLUI_Button		 *button_genFF, *button_doPR;
GLUI_Checkbox	 *cbox_showCurrentPatch, *cbox_showAmient, *cbox_smoothShade;
GLUI_Spinner	 *spinner_iterationLevel;

// IDs for callbacks
#define CB_UNSHOTPATCH_ID	100
#define CB_AMBIENT_ID		101
#define CB_SHADED_ID		102
#define SPIN_ITERATE_ID		103
#define BTN_GENFF		104
#define BTN_RUNPR		105

// Ambient term variables
Color reflectionFactor;	// overall interreflection factor R
Color ambient;		// total ambient factor
Color dAmbient;		// chance in ambience

// Global Control Variables
int mainWindow;
int totalStep 			= 0;
int currentPatchID;
int numOfIteration 		= 1;
int displayCurrentShotPatch 	= true;
int showAmbient 		= true;
int smoothShade 		= true;


//		Structures		//

// patch structure
typedef struct {
	int 	id;
	int 	vertices[4];	// array of vertices (index into VertexArray)
	Color 	emissivity;	// emitted radiosity
	Color 	reflectance;	// reflectance
	Vertex 	center;		// center of the patch
	Vector 	normal;		// patch normal
	double 	area;		// area of the patch
	Color 	radiosity;	// radiosity of the patch
	Color 	unshot;		// unshot radiosity of the patch
	int 	numelements;	// number of elements in patch
	int 	startelement;	// number of first element for this patch in ElementArray

} Patch;

// element structure
typedef struct {
	int 	id;
	int 	vertices[4];	// vertices of patch (index into VertexArray)
	Vertex 	center;		// center of the element
	double 	area;		// area of the element
	Patch* 	patch;		// Patch that this is an element of
	Color 	radiosity;	// radiosity of the element
} Element;

// unshot tag structure
// stores unshot, used for putting it in priority queue
struct UnshotTag {
	Color unshot;
	int id;

	UnshotTag(Color unshot_, int id_) : unshot(unshot_), id(id_) {}
} ;

// unshot tag comparison method for priority queue
struct CompareTag {
	// comparison operator
	// for priority queue
	bool operator ()(const UnshotTag &p1, const UnshotTag &p2) const {
		double p1Size = p1.unshot.r * p1.unshot.r + p1.unshot.g * p1.unshot.g + p1.unshot.b * p1.unshot.b;
		double p2Size = p2.unshot.r * p2.unshot.r + p2.unshot.g * p2.unshot.g + p2.unshot.b * p2.unshot.b;
		
		if (p1Size < p2Size){ return true; }
		else{ return false; }
	}
};

// Array buffers
int 	NumVertices;
Vertex* VertexArray;
Color* 	VertexColors;
int 	NumPatches;
Patch* 	PatchArray;
int 	NumElements;
Element* ElementArray;
double** lookUpTable;
priority_queue<UnshotTag, vector<UnshotTag>, CompareTag> unshotPatchQueue;

class Hemicube {

	struct HemiPixel {
		double formFactor;	// total form factor stored in pixel
		vector<int> id;		// list of element id that is projected

		HemiPixel() { formFactor = 0; }
	};

public: 
	vec3 center;
	vec3 normal;	// z-axis
	vec3 u, v;
	vector<HemiPixel*> hemiPixel;
	// construcors
	Hemicube() {};
	Hemicube(vec3 c, vec3 n, vec3 up) : center(c), normal(n), u(up) {

		// generate vector v
		v = cross(normal, u);
		// generate hemipixel set
		hemiPixel.push_back(new HemiPixel[HEMICUBE_SUBDIV * HEMICUBE_SUBDIV]);		// front face
		hemiPixel.push_back(new HemiPixel[HEMICUBE_SUBDIV * HEMICUBE_SUBDIV / 2]);	// side face
		hemiPixel.push_back(new HemiPixel[HEMICUBE_SUBDIV * HEMICUBE_SUBDIV / 2]);
		hemiPixel.push_back(new HemiPixel[HEMICUBE_SUBDIV * HEMICUBE_SUBDIV / 2]);
		hemiPixel.push_back(new HemiPixel[HEMICUBE_SUBDIV * HEMICUBE_SUBDIV / 2]);
	}
	~Hemicube() { for(int i=0; i<5; i++) delete[] hemiPixel[i]; }	// so important. or else it eats so much memory during runtime

	void computeDeltaFormFactor(Element element, int element_id) {
		for (int SIDE = 0; SIDE < 5; SIDE++) {

			int width, height;
			float left, right, bottom, top;
			vec3 lookat, up;

			// 1. set OpenGL viewport
			switch (SIDE) {
			case FRONT:
				width 	= height = HEMICUBE_SUBDIV;
				left 	= - HEMICUBE_HEIGHT;
				right 	= HEMICUBE_HEIGHT;
				bottom	= -HEMICUBE_HEIGHT;
				top 	= HEMICUBE_HEIGHT;
				lookat 	= center + normal;
				up 	= u;
				break;
			case TOP:
				width 	= HEMICUBE_SUBDIV;
				height 	= HEMICUBE_SUBDIV / 2;
				left 	= -HEMICUBE_HEIGHT;
				right 	= HEMICUBE_HEIGHT;
				bottom 	= 0;
				top 	= HEMICUBE_HEIGHT;
				lookat 	= center + u;
				up 	= normal;
				break;
			case RIGHT:
				width 	= HEMICUBE_SUBDIV;
				height 	= HEMICUBE_SUBDIV / 2;
				left 	= -HEMICUBE_HEIGHT;
				right 	= HEMICUBE_HEIGHT;
				bottom 	= 0;
				top 	= HEMICUBE_HEIGHT;
				lookat 	= center + v;
				up 	= normal;
				break;
			case BOTTOM:
				width 	= HEMICUBE_SUBDIV;
				height 	= HEMICUBE_SUBDIV / 2;
				left 	= -HEMICUBE_HEIGHT;
				right 	= HEMICUBE_HEIGHT;
				bottom 	= 0;
				top 	= HEMICUBE_HEIGHT;
				lookat 	= center - u;
				up 	= normal;
				break;
			case LEFT:
				width 	= HEMICUBE_SUBDIV;
				height 	= HEMICUBE_SUBDIV / 2;
				left 	= -HEMICUBE_HEIGHT;
				right 	= HEMICUBE_HEIGHT;
				bottom 	= 0;
				top 	= HEMICUBE_HEIGHT;
				lookat 	= center - v;
				up 	= normal;
				break;
			}

			// 2. render element to plane
			// element will be colored white

			glViewport(0, 0, width, height);
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glFrustum(left, right, bottom, top, HEMICUBE_HEIGHT, 10000);
			gluLookAt(center.x, center.y, center.z, lookat.x, lookat.y, lookat.z, up.x, up.y, up.z);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			// render current element white
			// everything else black.
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_QUADS);
			for (int id = 0; id < NumElements; id++) {
				
				if (id == element_id) glColor3f(1, 1, 1);
				else glColor3d(0, 0, 0);
				
				glVertex3f(VertexArray[ElementArray[id].vertices[0]].x, VertexArray[ElementArray[id].vertices[0]].y, VertexArray[ElementArray[id].vertices[0]].z);
				glVertex3f(VertexArray[ElementArray[id].vertices[1]].x, VertexArray[ElementArray[id].vertices[1]].y, VertexArray[ElementArray[id].vertices[1]].z);
				glVertex3f(VertexArray[ElementArray[id].vertices[2]].x, VertexArray[ElementArray[id].vertices[2]].y, VertexArray[ElementArray[id].vertices[2]].z);
				glVertex3f(VertexArray[ElementArray[id].vertices[3]].x, VertexArray[ElementArray[id].vertices[3]].y, VertexArray[ElementArray[id].vertices[3]].z);
			}
			glEnd();
			
			// read render data
			unsigned char* pixel = new unsigned char[width * height * 3];
			glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixel);

			// 3. find projected pixels
			//	  and save dFormFactor, update list of patch id in pixel
			for (int y = 0; y < height; y++) {			// for each row y
				for (int x = 0; x < width; x++) {		// each x in row
					int index = (y * width + x) * 3;

					if ((int)pixel[index] > 0) {	// if pixel is colored(projected)
						if (SIDE == FRONT) {
							// compute delta form factor
							double _y = (double)(y - height / 2) / HEMICUBE_SUBDIV;	// normalize to 1
							double _x = (double)(x - width / 2) / HEMICUBE_SUBDIV;	// normalize to 1
							double r = sqrt(_x*_x + _y*_y + HEMICUBE_HEIGHT);
							double dArea = 1.f / (HEMICUBE_SUBDIV * HEMICUBE_SUBDIV);
							double dFormFactor = (1 / (PI * pow(r , 4))) * dArea;	// THE LEGENDARY DELTA FORMFACTOR

							// add deltaFormfactor
							hemiPixel[SIDE][index/3].formFactor += dFormFactor;
							// add id
							hemiPixel[SIDE][index/3].id.push_back(element_id);
						}						
						else {
							// compute delta form factor
							double _y = (double)y / HEMICUBE_SUBDIV;				// normalize to 1
							double _x = (double)(x - width / 2) / HEMICUBE_SUBDIV;	// normalize to 1
							double r = sqrt(_x*_x + _y*_y + HEMICUBE_HEIGHT);
							double dArea = 1.f / (HEMICUBE_SUBDIV * HEMICUBE_SUBDIV /2);
							double dFormFactor = (_y / (PI * pow(r, 4))) * dArea;	// THE LEGENDARY DELTA FORMFACTOR

							// add deltaFormfactor
							hemiPixel[SIDE][index/3].formFactor += dFormFactor;
							// add id
							hemiPixel[SIDE][index/3].id.push_back(element_id);
						}
					}
				}
			}
			// delete memory.
			delete[] pixel;
		}
	}

	void updateLookUpTable(int patch_id) {
		
		// for every pixel,
		// pull id from vector
		// and 'add'(not assign) ff to table

		for (int SIDE = 0; SIDE < 5; SIDE++) {
			if (SIDE == FRONT) {
				for (int pix = 0; pix < HEMICUBE_SUBDIV*HEMICUBE_SUBDIV; pix++) {
					while (int size = hemiPixel[SIDE][pix].id.size() > 0) {
						int element_id = hemiPixel[SIDE][pix].id[size - 1];	// get last element (id)
						hemiPixel[SIDE][pix].id.pop_back();					// remove last element

						// update look up table
						lookUpTable[patch_id][element_id] += hemiPixel[SIDE][pix].formFactor;
					}
				}
			}
			else {
				for (int pix = 0; pix < HEMICUBE_SUBDIV*HEMICUBE_SUBDIV/2; pix++) {
					while (int size = hemiPixel[SIDE][pix].id.size() > 0) {
						int element_id = hemiPixel[SIDE][pix].id[size - 1];	// get last element (id)
						hemiPixel[SIDE][pix].id.pop_back();					// remove last element

						// update look up table
						lookUpTable[patch_id][element_id] += hemiPixel[SIDE][pix].formFactor;
					}
				}
			}
		}
	}
};



//		Functions		//

// update priority queue
void updatePriorityQueue() {

	// reset priority queue
	while (!unshotPatchQueue.empty())
		unshotPatchQueue.pop();

	for (int id = 0; id < NumPatches; id++) {
		UnshotTag unshot(PatchArray[id].unshot, id);
		unshotPatchQueue.push(unshot);
	}
}

// ** Do one step of Progressive Refinement ** //
void progressiveRefinement() {


	// get the most unshot patch (from priority queue)
	UnshotTag unshotTag = unshotPatchQueue.top();
	int mostUnshotID = unshotTag.id;
	currentPatchID = mostUnshotID;

	for (int element_id = 0; element_id < NumElements; element_id++) {

		// 1. determine increase in radiosity of element e due to Bi

		Color dRadiosity;

		Color reflectivity = ElementArray[element_id].patch->reflectance;
		Color unshot = PatchArray[mostUnshotID].unshot;
		double Fie = lookUpTable[mostUnshotID][element_id];
		double Ai = PatchArray[mostUnshotID].area;
		double Ae = ElementArray[element_id].area;
		double Aj = ElementArray[element_id].patch->area;

		dRadiosity.r = reflectivity.r * unshot.r * Fie * (Ai / Ae);
		dRadiosity.g = reflectivity.g * unshot.g * Fie * (Ai / Ae);
		dRadiosity.b = reflectivity.b * unshot.b * Fie * (Ai / Ae);


		// 2. add area weighted portion of of increased radiosity of element e
		//	  to radiosity of the patch j which contains element e

		Color *unshotJ = &ElementArray[element_id].patch->unshot;		// element's patch
		ElementArray[element_id].radiosity = ElementArray[element_id].radiosity + dRadiosity;	/// reflectivity * dAmbient; <- This(Ambient term) is just for display. Should not be added here
		unshotJ->r = unshotJ->r + dRadiosity.r*(Ae / Aj);
		unshotJ->g = unshotJ->g + dRadiosity.g*(Ae / Aj);
		unshotJ->b = unshotJ->b + dRadiosity.b*(Ae / Aj);

	}

	// 3. update Ambient

	double areaSum = 0;					// Sum of all patch area
	Color rfltWeightedAreaSum(0, 0, 0);	// Sum of all patch area * patch reflectance
	for (int p_id = 0; p_id < NumPatches; p_id++) {
		areaSum += PatchArray[p_id].area;

		rfltWeightedAreaSum.r += PatchArray[p_id].area * PatchArray[p_id].reflectance.r;
		rfltWeightedAreaSum.g += PatchArray[p_id].area * PatchArray[p_id].reflectance.g;
		rfltWeightedAreaSum.b += PatchArray[p_id].area * PatchArray[p_id].reflectance.b;
	}
	// overall interreflection factor R
	Color R, rfltAvg;
	rfltAvg.r = rfltWeightedAreaSum.r / areaSum;
	rfltAvg.g = rfltWeightedAreaSum.g / areaSum;
	rfltAvg.b = rfltWeightedAreaSum.b / areaSum;
	R.r = 1 / (1 - rfltAvg.r);
	R.g = 1 / (1 - rfltAvg.g);
	R.b = 1 / (1 - rfltAvg.b);
	// average area of unshot radiosity
	Color avgUnshot(0, 0, 0);
	for (int p_id = 0; p_id < NumPatches; p_id++) {
		avgUnshot.r += PatchArray[p_id].unshot.r * (PatchArray[p_id].area / areaSum);
		avgUnshot.g += PatchArray[p_id].unshot.g * (PatchArray[p_id].area / areaSum);
		avgUnshot.b += PatchArray[p_id].unshot.b * (PatchArray[p_id].area / areaSum);
	}

	dAmbient = R * avgUnshot;

	// 4. reset things
	PatchArray[mostUnshotID].unshot = Color(0, 0, 0);	// current patch's unshot <- 0
	updatePriorityQueue();								// update prioirty queue

														// print current step
	cout << "-----------------------------------------------------------------" << endl;
	cout << "\tPR::Current Step " << totalStep++ << endl;
	cout << "\tPR::Current Unshot Patch: " << PatchArray[mostUnshotID].id << endl;
	cout << "\tPR::new Ambient factor: " << dAmbient.r << ", " << dAmbient.g << ", " << dAmbient.b << endl;
}

// compute form factor of entire scene
void generateFormFactorTable() {

	// time variable
	time_t timer = time(0);

	cout << "\nGenFormFactors::Start generating patch - element form factors... " << endl;
	cout << "GenFormFactors::Start time: " << localtime(&timer)->tm_hour << ":" << localtime(&timer)->tm_min << ":" << localtime(&timer)->tm_sec << endl;

	// up vectors for hemicube
	vec3 up(0, 0, 1);		// towards upper direction of scene
	vec3 toCam(1, 0, 0);	// towards the camera
	Hemicube *hemicube;

	// for all patches
	for (int patch_id = 0; patch_id < NumPatches; patch_id++) {

		cout << "GenFormFactors::computing patch " << patch_id << "/" << NumPatches << "..." << endl;

		// create hemicube for the patch
		// depends on wall direction.
		int cosVal = dot(PatchArray[patch_id].normal, up);
		if (cosVal == 0)
			hemicube = new Hemicube(PatchArray[patch_id].center, PatchArray[patch_id].normal, up);
		else
			hemicube = new Hemicube(PatchArray[patch_id].center, PatchArray[patch_id].normal, toCam);

		for (int element_id = 0; element_id < NumElements; element_id++) {
			// compute form factor for that element
			hemicube->computeDeltaFormFactor(ElementArray[element_id], element_id);
		}
		// update look up table
		hemicube->updateLookUpTable(patch_id);

		// delete from memory
		delete hemicube;
	}

	// write file

	ofstream file;
	string fileName = "LookUpTable_output.csv";
	file.open(fileName);
	file << "F/E";
	for (int id = 0; id < NumElements; id++)	// write first row (element id)
		file << "," << id;
	file << endl;

	for (int p_id = 0; p_id < NumPatches; p_id++) {
		file << p_id;
		for (int e_id = 0; e_id < NumElements; e_id++) {
			file << "," << lookUpTable[p_id][e_id];
		}
		// move to next line
		file << endl;
	}
	file.close();

	timer = time(0);
	cout << "GenFormFactors::Done. Output file:" << fileName << endl;
	cout << "GenFormFactors::End time: " << localtime(&timer)->tm_hour << ":" << localtime(&timer)->tm_min << ":" << localtime(&timer)->tm_sec << endl;
}

// load file "LookUpTable.csv"
void loadLookUpTable() {

	cout << "LoadLUTable::start reading file..." << endl;

	ifstream file;
	double totalSum = 0;
	file.open("LookUpTable.csv");

	// read first line - element id
	string firstLine;
	getline(file, firstLine);
	// read table
	for (int p_id = 0; p_id < NumPatches; p_id++) {

		// read first column of the row (patch id)
		int id;
		file >> id;
		for (int e_id = 0; e_id < NumElements; e_id++) {
			char comma;
			file >> comma;

			double val;
			file >> val;
			lookUpTable[p_id][e_id] = val;
			totalSum += val;
		}
	}
	file.close();

	cout << "LoadLUTable::total form factor: " << totalSum << endl;
	cout << "LoadLUTable::file reading done.\n" << endl;
}

// update vertex color, add amience
void updateVertexColor() {

	// bilinear interpolation
	// with neighboring patches

	struct ColorStack {
		int count;
		Color color;
		ColorStack() : count(0), color(Color(0, 0, 0)) {}
	};

	ColorStack *colorStack = new ColorStack[NumVertices];

	// update vertex radiosity
	// by each element
	for (int i = 0; i < NumElements; i++) {

		// add radiosity color

		if (showAmbient) {
			colorStack[ElementArray[i].vertices[0]].color += ElementArray[i].radiosity + ElementArray[i].patch->reflectance * dAmbient;
			colorStack[ElementArray[i].vertices[1]].color += ElementArray[i].radiosity + ElementArray[i].patch->reflectance * dAmbient;
			colorStack[ElementArray[i].vertices[2]].color += ElementArray[i].radiosity + ElementArray[i].patch->reflectance * dAmbient;
			colorStack[ElementArray[i].vertices[3]].color += ElementArray[i].radiosity + ElementArray[i].patch->reflectance * dAmbient;
		}
		else {
			colorStack[ElementArray[i].vertices[0]].color += ElementArray[i].radiosity;
			colorStack[ElementArray[i].vertices[1]].color += ElementArray[i].radiosity;
			colorStack[ElementArray[i].vertices[2]].color += ElementArray[i].radiosity;
			colorStack[ElementArray[i].vertices[3]].color += ElementArray[i].radiosity;
		}
		// increment count
		colorStack[ElementArray[i].vertices[0]].count += 1;
		colorStack[ElementArray[i].vertices[1]].count += 1;
		colorStack[ElementArray[i].vertices[2]].count += 1;
		colorStack[ElementArray[i].vertices[3]].count += 1;
	}

	// average color
	for (int i = 0; i < NumVertices; i++) {
		VertexColors[i].r = colorStack[i].color.r / colorStack[i].count;
		VertexColors[i].g = colorStack[i].color.g / colorStack[i].count;
		VertexColors[i].b = colorStack[i].color.b / colorStack[i].count;
	}

	delete colorStack;
}

// draw patch by id
void drawPatch(int id) {
	glColor3f(0, 1, 1);
	glPointSize(3.f);
	glBegin(GL_LINE_STRIP);
	glVertex3f(VertexArray[PatchArray[id].vertices[0]].x, VertexArray[PatchArray[id].vertices[0]].y, VertexArray[PatchArray[id].vertices[0]].z);
	glVertex3f(VertexArray[PatchArray[id].vertices[1]].x, VertexArray[PatchArray[id].vertices[1]].y, VertexArray[PatchArray[id].vertices[1]].z);
	glVertex3f(VertexArray[PatchArray[id].vertices[2]].x, VertexArray[PatchArray[id].vertices[2]].y, VertexArray[PatchArray[id].vertices[2]].z);
	glVertex3f(VertexArray[PatchArray[id].vertices[3]].x, VertexArray[PatchArray[id].vertices[3]].y, VertexArray[PatchArray[id].vertices[3]].z);
	glVertex3f(VertexArray[PatchArray[id].vertices[0]].x, VertexArray[PatchArray[id].vertices[0]].y, VertexArray[PatchArray[id].vertices[0]].z);
	glEnd();
}

// draw element by id
void drawElement(int id) {
	glColor3f(1, 0.5, 0.5);
	glPointSize(2.f);
	glBegin(GL_QUADS);
	glVertex3f(VertexArray[ElementArray[id].vertices[0]].x, VertexArray[ElementArray[id].vertices[0]].y, VertexArray[ElementArray[id].vertices[0]].z);
	glVertex3f(VertexArray[ElementArray[id].vertices[1]].x, VertexArray[ElementArray[id].vertices[1]].y, VertexArray[ElementArray[id].vertices[1]].z);
	glVertex3f(VertexArray[ElementArray[id].vertices[2]].x, VertexArray[ElementArray[id].vertices[2]].y, VertexArray[ElementArray[id].vertices[2]].z);
	glVertex3f(VertexArray[ElementArray[id].vertices[3]].x, VertexArray[ElementArray[id].vertices[3]].y, VertexArray[ElementArray[id].vertices[3]].z);
	glVertex3f(VertexArray[ElementArray[id].vertices[0]].x, VertexArray[ElementArray[id].vertices[0]].y, VertexArray[ElementArray[id].vertices[0]].z);
	glEnd();
}

// initialize scene : compute initial ambience, init unshot patches
void initScene(){

	// 0. compute some factors...

	Color avgPatchRefl(0, 0, 0); // weighted average of the patch reflectives
	Color sumEmi(0, 0, 0);		 // weighted sum of emission

	double areaSum = 0;
	for (int id = 0; id < NumPatches; id++) {

		avgPatchRefl.r += (PatchArray[id].reflectance.r * PatchArray[id].area);
		avgPatchRefl.g += (PatchArray[id].reflectance.g * PatchArray[id].area);
		avgPatchRefl.b += (PatchArray[id].reflectance.b * PatchArray[id].area);

		sumEmi.r += (PatchArray[id].emissivity.r * PatchArray[id].area);
		sumEmi.g += (PatchArray[id].emissivity.g * PatchArray[id].area);
		sumEmi.b += (PatchArray[id].emissivity.b * PatchArray[id].area);

		areaSum += PatchArray[id].area;
	}


	// 1. compute reflection factor R

	avgPatchRefl.r /= areaSum;
	avgPatchRefl.g /= areaSum;
	avgPatchRefl.b /= areaSum;
	// assign final value
	reflectionFactor.r = 1 / (1 - avgPatchRefl.r);
	reflectionFactor.g = 1 / (1 - avgPatchRefl.g);
	reflectionFactor.b = 1 / (1 - avgPatchRefl.b);

	cout << "Init::reflection factor R: " << reflectionFactor.r << ", " << reflectionFactor.g << ", " << reflectionFactor.b << endl;


	// 2. determine initial ambient from given emission

	ambient.r = reflectionFactor.r * sumEmi.r / areaSum;
	ambient.g = reflectionFactor.g * sumEmi.g / areaSum;
	ambient.b = reflectionFactor.b * sumEmi.b / areaSum;

	dAmbient.r = 0;
	dAmbient.g = 0;
	dAmbient.b = 0;

	cout << "Init::ambient factor: " << ambient.r << ", " << ambient.g << ", " << ambient.b << endl;

	// 3. initialize unshot radiosity with emission value

	for (int id = 0; id < NumPatches; id++) {
		PatchArray[id].unshot = PatchArray[id].emissivity;
	}

	// 4. initialize look up table

	// patch to element table
	lookUpTable = new double*[NumPatches];
	for (int p_id = 0; p_id < NumPatches; p_id++) {
		lookUpTable[p_id] = new double[NumElements];
		for (int e_id = 0; e_id < NumElements; e_id++)
			lookUpTable[p_id][e_id] = 0;		// initialize value
	}

	// 5. update initial heap
	//	  & initial vertex color
	updateVertexColor();
	updatePriorityQueue();
}

int loadData(void) {

	int i, j, k;
	int nverts, vertnum, startvert;
	int elnum;
	Vertex* vtemp;
	std::ifstream infi("scene.dat");
	Vector v1;
	Vector v2;
	double length1, length2;
	double temparea;

	NumElements = 0;

	// read initial vertices
	infi >> nverts;
	vtemp = new Vertex[nverts];
	for (i = 0; i<nverts; i++) {
		infi >> vtemp[i].x >> vtemp[i].y >> vtemp[i].z;
	}
	NumVertices = nverts;

	// read patches
	infi >> NumPatches;
	PatchArray = new Patch[NumPatches];
	for (i = 0; i<NumPatches; i++) {

		// Patch id
		PatchArray[i].id = i;

		// Read patch i
		infi >> PatchArray[i].vertices[0] >> PatchArray[i].vertices[1]
			>> PatchArray[i].vertices[2] >> PatchArray[i].vertices[3];
		infi >> PatchArray[i].emissivity.r >> PatchArray[i].emissivity.g
			>> PatchArray[i].emissivity.b;
		infi >> PatchArray[i].reflectance.r >> PatchArray[i].reflectance.g
			>> PatchArray[i].reflectance.b;
		infi >> PatchArray[i].numelements;

		// **** TEST : divide more **** //
		PatchArray[i].numelements += 2;

		if (PatchArray[i].emissivity.r > 0 || PatchArray[i].emissivity.g > 0 || PatchArray[i].emissivity.b > 0){
			cout << "Load::incident light patch " << i << endl;
		}

		NumVertices += (PatchArray[i].numelements + 1) *
			(PatchArray[i].numelements + 1);
		PatchArray[i].startelement = NumElements;
		NumElements += (PatchArray[i].numelements * PatchArray[i].numelements);

		// patch center
		PatchArray[i].center.x = (vtemp[PatchArray[i].vertices[0]].x +
			vtemp[PatchArray[i].vertices[1]].x +
			vtemp[PatchArray[i].vertices[2]].x +
			vtemp[PatchArray[i].vertices[3]].x) / 4.0;
		PatchArray[i].center.y = (vtemp[PatchArray[i].vertices[0]].y +
			vtemp[PatchArray[i].vertices[1]].y +
			vtemp[PatchArray[i].vertices[2]].y +
			vtemp[PatchArray[i].vertices[3]].y) / 4.0;
		PatchArray[i].center.z = (vtemp[PatchArray[i].vertices[0]].z +
			vtemp[PatchArray[i].vertices[1]].z +
			vtemp[PatchArray[i].vertices[2]].z +
			vtemp[PatchArray[i].vertices[3]].z) / 4.0;

		// patch area
		v1.x = vtemp[PatchArray[i].vertices[1]].x -
			vtemp[PatchArray[i].vertices[0]].x;
		v1.y = vtemp[PatchArray[i].vertices[1]].y -
			vtemp[PatchArray[i].vertices[0]].y;
		v1.z = vtemp[PatchArray[i].vertices[1]].z -
			vtemp[PatchArray[i].vertices[0]].z;
		v2.x = vtemp[PatchArray[i].vertices[3]].x -
			vtemp[PatchArray[i].vertices[0]].x;
		v2.y = vtemp[PatchArray[i].vertices[3]].y -
			vtemp[PatchArray[i].vertices[0]].y;
		v2.z = vtemp[PatchArray[i].vertices[3]].z -
			vtemp[PatchArray[i].vertices[0]].z;
		length1 = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
		length2 = sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
		PatchArray[i].area = length1*length2;
		v1.x /= length1;
		v1.y /= length1;
		v1.z /= length1;
		v2.x /= length2;
		v2.y /= length2;
		v2.z /= length2;
		PatchArray[i].normal.x = v1.y*v2.z - v1.z*v2.y;
		PatchArray[i].normal.y = v2.x*v1.z - v1.x*v2.z;
		PatchArray[i].normal.z = v1.x*v2.y - v1.y*v2.x;
		PatchArray[i].radiosity.r = PatchArray[i].emissivity.r;
		PatchArray[i].radiosity.g = PatchArray[i].emissivity.g;
		PatchArray[i].radiosity.b = PatchArray[i].emissivity.b;
		PatchArray[i].unshot.r = PatchArray[i].emissivity.r;
		PatchArray[i].unshot.g = PatchArray[i].emissivity.g;
		PatchArray[i].unshot.b = PatchArray[i].emissivity.b;
	}

	// create elements (including new vertices)
	VertexArray = new Vertex[NumVertices];
	VertexColors = new Color[NumVertices];
	ElementArray = new Element[NumElements];

	// Copy original (patch) vertices to beginning of array
	for (i = 0; i<nverts; i++) {
		VertexArray[i].x = vtemp[i].x;
		VertexArray[i].y = vtemp[i].y;
		VertexArray[i].z = vtemp[i].z;
		VertexColors[i].r = 0.0;
		VertexColors[i].g = 0.0;
		VertexColors[i].b = 0.0;
	}

	// Form Vertices for new elements
	vertnum = nverts;
	elnum = 0;
	for (i = 0; i<NumPatches; i++) {
		v1.x = (VertexArray[PatchArray[i].vertices[1]].x -
			VertexArray[PatchArray[i].vertices[0]].x) /
			PatchArray[i].numelements;
		v1.y = (VertexArray[PatchArray[i].vertices[1]].y -
			VertexArray[PatchArray[i].vertices[0]].y) /
			PatchArray[i].numelements;
		v1.z = (VertexArray[PatchArray[i].vertices[1]].z -
			VertexArray[PatchArray[i].vertices[0]].z) /
			PatchArray[i].numelements;
		v2.x = (VertexArray[PatchArray[i].vertices[3]].x -
			VertexArray[PatchArray[i].vertices[0]].x) /
			PatchArray[i].numelements;
		v2.y = (VertexArray[PatchArray[i].vertices[3]].y -
			VertexArray[PatchArray[i].vertices[0]].y) /
			PatchArray[i].numelements;
		v2.z = (VertexArray[PatchArray[i].vertices[3]].z -
			VertexArray[PatchArray[i].vertices[0]].z) /
			PatchArray[i].numelements;

		startvert = vertnum;

		for (j = 0; j<PatchArray[i].numelements + 1; j++) {
			for (k = 0; k<PatchArray[i].numelements + 1; k++) {
				// Create new vertex
				VertexArray[vertnum].x =
					VertexArray[PatchArray[i].vertices[0]].x +
					k * v1.x + j * v2.x;
				VertexArray[vertnum].y =
					VertexArray[PatchArray[i].vertices[0]].y +
					k * v1.y + j * v2.y;
				VertexArray[vertnum].z =
					VertexArray[PatchArray[i].vertices[0]].z +
					k * v1.z + j * v2.z;
				VertexColors[vertnum].r = 0.0;
				VertexColors[vertnum].g = 0.0;
				VertexColors[vertnum].b = 0.0;

				vertnum++;
			}
		}
		temparea = PatchArray[i].area /
			(PatchArray[i].numelements * PatchArray[i].numelements);
		// Form Elements for new patch
		for (j = 0; j<PatchArray[i].numelements; j++) {
			for (k = 0; k<PatchArray[i].numelements; k++) {

				// Element id
				ElementArray[elnum].id = elnum;

				// Set vertices
				ElementArray[elnum].vertices[0] = startvert + k +
					j * (PatchArray[i].numelements + 1);
				ElementArray[elnum].vertices[1] = startvert + (k + 1) +
					j * (PatchArray[i].numelements + 1);
				ElementArray[elnum].vertices[2] = startvert + (k + 1) +
					(j + 1) * (PatchArray[i].numelements + 1);
				ElementArray[elnum].vertices[3] = startvert + k +
					(j + 1) * (PatchArray[i].numelements + 1);

				// Set center
				ElementArray[elnum].center.x =
					(VertexArray[ElementArray[elnum].vertices[0]].x +
						VertexArray[ElementArray[elnum].vertices[1]].x +
						VertexArray[ElementArray[elnum].vertices[2]].x +
						VertexArray[ElementArray[elnum].vertices[3]].x) / 4.0;
				ElementArray[elnum].center.y =
					(VertexArray[ElementArray[elnum].vertices[0]].y +
						VertexArray[ElementArray[elnum].vertices[1]].y +
						VertexArray[ElementArray[elnum].vertices[2]].y +
						VertexArray[ElementArray[elnum].vertices[3]].y) / 4.0;
				ElementArray[elnum].center.z =
					(VertexArray[ElementArray[elnum].vertices[0]].z +
						VertexArray[ElementArray[elnum].vertices[1]].z +
						VertexArray[ElementArray[elnum].vertices[2]].z +
						VertexArray[ElementArray[elnum].vertices[3]].z) / 4.0;

				ElementArray[elnum].area = temparea;
				ElementArray[elnum].patch = &PatchArray[i];
				ElementArray[elnum].radiosity = PatchArray[i].radiosity;
				elnum++;
			}
		}
	}

	delete[] vtemp;
	infi.close();
	return 0;
}

// initialize entire scene
void init(void) {
	cout << "----------------------------------------" << endl;
	cout << "\n\tInitializeing Scene...\n" << endl;

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnable(GL_DEPTH_TEST);
	loadData();		// load scene data
	initScene();	// init initial scene factors
	loadLookUpTable();

	cout << "\n\tInitialization Complete\n" << endl;
	cout << "----------------------------------------" << endl;
}

// GLUT display
void display(void) {
	int i;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (displayCurrentShotPatch)
		drawPatch(currentPatchID);

	glEnable(GL_FLAT | GL_SMOOTH);
	if (smoothShade)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	glBegin(GL_QUADS);
	for (i = 0; i<NumElements; i++) {

		glColor3f(VertexColors[ElementArray[i].vertices[0]].r,
			VertexColors[ElementArray[i].vertices[0]].g,
			VertexColors[ElementArray[i].vertices[0]].b);
		glVertex3f(VertexArray[ElementArray[i].vertices[0]].x,
			VertexArray[ElementArray[i].vertices[0]].y,
			VertexArray[ElementArray[i].vertices[0]].z);
		glColor3f(VertexColors[ElementArray[i].vertices[1]].r,
			VertexColors[ElementArray[i].vertices[1]].g,
			VertexColors[ElementArray[i].vertices[1]].b);
		glVertex3f(VertexArray[ElementArray[i].vertices[1]].x,
			VertexArray[ElementArray[i].vertices[1]].y,
			VertexArray[ElementArray[i].vertices[1]].z);
		glColor3f(VertexColors[ElementArray[i].vertices[2]].r,
			VertexColors[ElementArray[i].vertices[2]].g,
			VertexColors[ElementArray[i].vertices[2]].b);
		glVertex3f(VertexArray[ElementArray[i].vertices[2]].x,
			VertexArray[ElementArray[i].vertices[2]].y,
			VertexArray[ElementArray[i].vertices[2]].z);
		glColor3f(VertexColors[ElementArray[i].vertices[3]].r,
			VertexColors[ElementArray[i].vertices[3]].g,
			VertexColors[ElementArray[i].vertices[3]].b);
		glVertex3f(VertexArray[ElementArray[i].vertices[3]].x,
			VertexArray[ElementArray[i].vertices[3]].y,
			VertexArray[ElementArray[i].vertices[3]].z);
	}
	glEnd();
	glFlush();
}

// GLUT reshape
void reshape(int w, int h) {
	
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-5.0, 5.0, -4.0, 4.0, 10.1, 25.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(20.0, 5.0, 4.0,
		10.0, 5.0, 4.0,
		0.0, 0.0, 1.0);
}

void buttonCallback(GLUI_Control* control) {

	if (control->get_id() == BTN_RUNPR) {
		for (int i = 0; i < numOfIteration; i++) {
			progressiveRefinement();
		}
		updateVertexColor();
		glutPostRedisplay();
	}
	else if (control->get_id() == BTN_GENFF) {
		generateFormFactorTable();		// compute form factors
	}
}

int main(int argc, char** argv)
{
	cout << "Project3 - Computer Graphics, Fall 2017, Texas A&M University" << endl;
	cout << "Made by Somyung (David) Oh.\n" << endl;

	// GLUT initialization
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(800, 600);
	glutInitWindowPosition(100, 100);
	mainWindow = glutCreateWindow("Radiosity Solver");
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);

	// GLUI initialization
	glui = GLUI_Master.create_glui("Render Control", 0, 800, 70);
	panel_control 			= new  GLUI_Panel(glui, "Settings");
	cbox_showCurrentPatch 	= new GLUI_Checkbox(panel_control, "show current patch", &displayCurrentShotPatch, CB_UNSHOTPATCH_ID, buttonCallback);
	cbox_showAmient 		= new GLUI_Checkbox(panel_control, "show ambient light", &showAmbient,CB_AMBIENT_ID, buttonCallback);
	cbox_smoothShade 		= new GLUI_Checkbox(panel_control, "gouraud shade", &smoothShade, CB_SHADED_ID, buttonCallback);
	spinner_iterationLevel 	= new GLUI_Spinner(panel_control, "interation in step", &numOfIteration, -1, buttonCallback);
	spinner_iterationLevel->set_int_limits(1, 150, GLUI_LIMIT_CLAMP);
	spinner_iterationLevel->set_speed(0.05);
	button_doPR 			= new GLUI_Button(glui, "Do Progressive Refinement", BTN_RUNPR, buttonCallback);
	glui->add_separator();
	button_genFF 			= new GLUI_Button(glui, "Generate Form Factor", BTN_GENFF, buttonCallback);
	glui->set_main_gfx_window(mainWindow);

	glutMainLoop();

	return 0;
}
