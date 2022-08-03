#ifndef OBJECT_H
#define OBJECT_H
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <time.h>
#include "vector.h"
#define BUFFER_LENGTH 256
#define LARGEST_DISTANCE 512
#define LENGTH 800
#define HEIGHT 600
#define VIEWPORT_DISTANCE 0.01
#define FACET 3
#define DIM 3
#define ER 1.0e-3

class Vertex {
public:
	std::vector<double> coordinate{ std::vector<double>(3) };
public:
	Vertex() {}
	Vertex(const double x, const double y, const double z = 0) {
		coordinate[0] = x;
		coordinate[1] = y;
		coordinate[2] = z;
	}
	~Vertex() {}
	Vertex(const Vertex& v) {
		std::vector<double> vCoord(3);
		v.getCoordinate(&vCoord);
		if (!setCoordinate(vCoord))
			std::cout << "ERROR::Invalid Vertex Copy Constructing!" << std::endl;
	}
	int getCoordinate(std::vector<double>* coord) const {
		if ((*coord).size() < 3)
			(*coord).resize(3);
		(*coord)[0] = coordinate[0];
		(*coord)[1] = coordinate[1];
		(*coord)[2] = coordinate[2];
		return 1;
	}
	int setCoordinate(std::vector<double> coord) {
		if (coord.size() < 3)
			return -1;
		coordinate[0] = coord[0];
		coordinate[1] = coord[1];
		coordinate[2] = coord[2];
		return 1;
	}
	double x() const {
		return coordinate[0];
	}
	double y() const {
		return coordinate[1];
	}
	double z() const {
		return coordinate[2];
	}
};

class Object {
protected:
	std::vector<double> vertices;
	std::vector< std::vector<int> > polygons;
	std::set< std::set<int> > edges;
	std::vector<std::vector<int> > edgesV;
	std::vector<double> newVertices;
	std::vector< std::vector<int> > newPolygons;

	std::vector<unsigned short int> colors;
	std::vector<int> triangles;
	std::vector<unsigned short> triColors;
	std::vector<double> verticesUVNs;
	std::vector<double> verticesUPVPs;
	int establishVerticesUVNs(Vertex viewPoint, Vector viewDirection);
	int establishVerticesUPVPs(double distance);
	int clearFrameBuffer();
	int clearZBuffer();
public:
	unsigned short frameBuffer[LENGTH*HEIGHT * 3];
	double zBuffer[LENGTH*HEIGHT];
	int ZBuffer(Vertex viewPoint, Vector viewDirection);
	int writePpmFile(const char* filename);
	
public:
	Object() {};
	~Object() {};
	int readObjFile(const char* filename);
	int writeObjFile(const char* filename);
	int DooSabin(int times = 1);
	int CatmullClark(int times = 1);
	int Loop(int times = 1);
	int establishEdges();
	int establishPolygonOffset();
	int establishVertex2Polygons();
	int establishEdge2Polygons();
	int establishVertex2Edges();
	int establishEdgesV();
	int establishPolygon2Edges();
	
protected:
	std::vector< std::vector<int> > vertex2Polygons;
	std::vector< std::vector<int> > edge2Polygons;
	std::vector< std::vector<int> > vertex2Edges;
	std::vector< std::set<int> > polygon2Edges;
	std::vector<int> polygonOffset;
public:
	int setColors();
	int tri();
};

#endif