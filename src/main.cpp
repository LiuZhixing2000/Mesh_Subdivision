//#define DOO_SABIN
//#define CATMULL_CLARK
#define LOOP
#define ZBUFFER
#include <iostream>
#include "object.h"
#include "glut.h"
short int colors[HEIGHT*LENGTH * 3];
void paint(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1, 1, 1, 1);
	glViewport(0, 0, LENGTH, HEIGHT);
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < LENGTH; j++) {
			float xx = (float)(j - LENGTH / 2) / (float)(LENGTH / 2);
			float yy = (float)(-i + HEIGHT / 2 - 1) / (float)(HEIGHT / 2);
			float rr = (float)colors[(LENGTH*i + j) * 3] / 255.0;
			float gg = (float)colors[(LENGTH*i + j) * 3 + 1] / 255.0;
			float bb = (float)colors[(LENGTH*i + j) * 3 + 2] / 255.0;
			glColor3f(rr, gg, bb);
			glBegin(GL_POINTS);
			glVertex2f(xx, yy);
			glEnd();
		}
	}
	glFlush();
	std::cout << "Visible Surface Determination Result visualization Success!" << std::endl;
	std::cout << "Click close button to end this execution." << std::endl;
}

int main(int argc, char** argv) {
	if (argc < 3) {
		std::cout << "ERROR::Invalid param in!" << std::endl;
		return -1;
	}
	std::string inFilename = std::string(argv[1]);
	std::string outFilename = inFilename;
	outFilename = outFilename.substr(0, outFilename.size() - 4);
	int times = argv[2][0] - '0';
	// 读取模型
	Object* obj = new Object;
	obj->readObjFile(inFilename.c_str());
	// 进行曲面细分
#ifdef DOO_SABIN
	obj->DooSabin(times);
	outFilename.append("_DooSabin_");
#endif
#ifdef CATMULL_CLARK
	obj->CatmullClark(times);
	outFilename.append("_CatmullClark_");
#endif
#ifdef LOOP
	obj->Loop(times);
	outFilename.append("_Loop_");
#endif
	
	// 写到obj文件
	char t[2] = { times + '0', '\0' };
	outFilename.append(std::string(t));
	outFilename.append(".obj");
	obj->writeObjFile(outFilename.c_str());

	//std::string commandLine = "..\\ParaView\\bin\\paraview.exe ";
	//commandLine.append(outFilename);
	//system(commandLine.c_str());
	// visualization
	obj->setColors();
	obj->tri();
	Vertex viewPoint;
	viewPoint.setCoordinate(std::vector<double>({ 0, 0, 5 }));
	Vector viewDirection(0, 0, -1);
	obj->ZBuffer(viewPoint, viewDirection);
	std::string ppmFilename = outFilename;
	ppmFilename = ppmFilename.substr(0, ppmFilename.size() - 4);
	ppmFilename.append("_ZBuffer.ppm");
	obj->writePpmFile(ppmFilename.c_str());
	std::cout << "Visible Surface Determination Result visualizing, please wait..." << std::endl;
	for (int i = 0; i < HEIGHT*LENGTH * 3; i++) {
		colors[i] = obj->frameBuffer[i];
	}
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowSize(LENGTH, HEIGHT);
	glutInitWindowPosition(500, 300);
#ifdef ZBUFFER
	glutCreateWindow("ZBuffer Visible Surface Determination Result");
#endif
#ifdef SCANLINE_ZBUFFER
	glutCreateWindow("ScanLineZBuffer Visible Surface Determination Result");
#endif
	glutDisplayFunc(paint);
	glutMainLoop();
	return 0;
}