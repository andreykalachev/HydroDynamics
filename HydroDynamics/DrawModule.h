#ifndef DRAWMODULE_H
#define DRAWMODULE_H

#include "GL/glut.h"
#include "fade3d/include_fade3d/Fade_3D.h"

using namespace std;
using namespace FADE3D;

vector<Tet3*> tets;

void displayCorners()
{
	glBegin(GL_LINES);

	// draws the edges of the tetrahedron
	for (int i = 0; i < tets.size(); i++)
	{

		glVertex3d(tets[i]->getCorner(0)->x(), tets[i]->getCorner(0)->y(), tets[i]->getCorner(0)->z());
		glVertex3d(tets[i]->getCorner(1)->x(), tets[i]->getCorner(1)->y(), tets[i]->getCorner(1)->z());

		glVertex3d(tets[i]->getCorner(0)->x(), tets[i]->getCorner(0)->y(), tets[i]->getCorner(0)->z());
		glVertex3d(tets[i]->getCorner(2)->x(), tets[i]->getCorner(2)->y(), tets[i]->getCorner(2)->z());

		glVertex3d(tets[i]->getCorner(0)->x(), tets[i]->getCorner(0)->y(), tets[i]->getCorner(0)->z());
		glVertex3d(tets[i]->getCorner(3)->x(), tets[i]->getCorner(3)->y(), tets[i]->getCorner(3)->z());

		glVertex3d(tets[i]->getCorner(1)->x(), tets[i]->getCorner(1)->y(), tets[i]->getCorner(1)->z());
		glVertex3d(tets[i]->getCorner(2)->x(), tets[i]->getCorner(2)->y(), tets[i]->getCorner(2)->z());

		glVertex3d(tets[i]->getCorner(1)->x(), tets[i]->getCorner(1)->y(), tets[i]->getCorner(1)->z());
		glVertex3d(tets[i]->getCorner(3)->x(), tets[i]->getCorner(3)->y(), tets[i]->getCorner(3)->z());

		glVertex3d(tets[i]->getCorner(2)->x(), tets[i]->getCorner(2)->y(), tets[i]->getCorner(2)->z());
		glVertex3d(tets[i]->getCorner(3)->x(), tets[i]->getCorner(3)->y(), tets[i]->getCorner(3)->z());

	}

	glEnd();

}

void displayVectorsB()
{

	glBegin(GL_LINES);

	for (int i = 0; i < tets.size(); i++)
	{
		double x, y, z, x_m, y_m, z_m;

		//b1

		x = (tets[i]->getCorner(2)->y() - tets[i]->getCorner(1)->y())*(tets[i]->getCorner(3)->z() - tets[i]->getCorner(1)->z()) -
			(tets[i]->getCorner(3)->y() - tets[i]->getCorner(1)->y())*(tets[i]->getCorner(2)->z() - tets[i]->getCorner(1)->z());

		y = (tets[i]->getCorner(2)->z() - tets[i]->getCorner(1)->z())*(tets[i]->getCorner(3)->x() - tets[i]->getCorner(1)->x()) -
			(tets[i]->getCorner(3)->z() - tets[i]->getCorner(1)->z())*(tets[i]->getCorner(2)->x() - tets[i]->getCorner(1)->x());

		z = (tets[i]->getCorner(2)->x() - tets[i]->getCorner(1)->x())*(tets[i]->getCorner(3)->y() - tets[i]->getCorner(1)->y()) -
			(tets[i]->getCorner(3)->x() - tets[i]->getCorner(1)->x())*(tets[i]->getCorner(2)->y() - tets[i]->getCorner(1)->y());

		x_m = (tets[i]->getCorner(1)->x() + tets[i]->getCorner(2)->x() + tets[i]->getCorner(3)->x()) / 3;
		y_m = (tets[i]->getCorner(1)->y() + tets[i]->getCorner(2)->y() + tets[i]->getCorner(3)->y()) / 3;
		z_m = (tets[i]->getCorner(1)->z() + tets[i]->getCorner(2)->z() + tets[i]->getCorner(3)->z()) / 3;

		glVertex3d(x_m, y_m, z_m);
		glVertex3d(x_m + x / 200, y_m + y / 200, z_m + z / 200);

		//b2


		x = (tets[i]->getCorner(3)->y() - tets[i]->getCorner(0)->y())*(tets[i]->getCorner(2)->z() - tets[i]->getCorner(0)->z()) -
			(tets[i]->getCorner(2)->y() - tets[i]->getCorner(0)->y())*(tets[i]->getCorner(3)->z() - tets[i]->getCorner(0)->z());

		y = (tets[i]->getCorner(3)->z() - tets[i]->getCorner(0)->z())*(tets[i]->getCorner(2)->x() - tets[i]->getCorner(0)->x()) -
			(tets[i]->getCorner(2)->z() - tets[i]->getCorner(0)->z())*(tets[i]->getCorner(3)->x() - tets[i]->getCorner(0)->x());

		z = (tets[i]->getCorner(3)->x() - tets[i]->getCorner(0)->x())*(tets[i]->getCorner(2)->y() - tets[i]->getCorner(0)->y()) -
			(tets[i]->getCorner(2)->x() - tets[i]->getCorner(0)->x())*(tets[i]->getCorner(3)->y() - tets[i]->getCorner(0)->y());

		x_m = (tets[i]->getCorner(0)->x() + tets[i]->getCorner(2)->x() + tets[i]->getCorner(3)->x()) / 3;
		y_m = (tets[i]->getCorner(0)->y() + tets[i]->getCorner(2)->y() + tets[i]->getCorner(3)->y()) / 3;
		z_m = (tets[i]->getCorner(0)->z() + tets[i]->getCorner(2)->z() + tets[i]->getCorner(3)->z()) / 3;

		glVertex3d(x_m, y_m, z_m);
		glVertex3d(x_m + x / 200, y_m + y / 200, z_m + z / 200);

		//b3

		x = (tets[i]->getCorner(0)->y() - tets[i]->getCorner(3)->y())*(tets[i]->getCorner(1)->z() - tets[i]->getCorner(3)->z()) -
			(tets[i]->getCorner(1)->y() - tets[i]->getCorner(3)->y())*(tets[i]->getCorner(0)->z() - tets[i]->getCorner(3)->z());

		y = (tets[i]->getCorner(0)->z() - tets[i]->getCorner(3)->z())*(tets[i]->getCorner(1)->x() - tets[i]->getCorner(3)->x()) -
			(tets[i]->getCorner(1)->z() - tets[i]->getCorner(3)->z())*(tets[i]->getCorner(0)->x() - tets[i]->getCorner(3)->x());

		z = (tets[i]->getCorner(0)->x() - tets[i]->getCorner(3)->x())*(tets[i]->getCorner(1)->y() - tets[i]->getCorner(3)->y()) -
			(tets[i]->getCorner(1)->x() - tets[i]->getCorner(3)->x())*(tets[i]->getCorner(0)->y() - tets[i]->getCorner(3)->y());

		x_m = (tets[i]->getCorner(0)->x() + tets[i]->getCorner(1)->x() + tets[i]->getCorner(3)->x()) / 3;
		y_m = (tets[i]->getCorner(0)->y() + tets[i]->getCorner(1)->y() + tets[i]->getCorner(3)->y()) / 3;
		z_m = (tets[i]->getCorner(0)->z() + tets[i]->getCorner(1)->z() + tets[i]->getCorner(3)->z()) / 3;

		glVertex3d(x_m, y_m, z_m);
		glVertex3d(x_m + x / 200, y_m + y / 200, z_m + z / 200);

		//b4

		x = (tets[i]->getCorner(1)->y() - tets[i]->getCorner(2)->y())*(tets[i]->getCorner(0)->z() - tets[i]->getCorner(2)->z()) -
			(tets[i]->getCorner(0)->y() - tets[i]->getCorner(2)->y())*(tets[i]->getCorner(1)->z() - tets[i]->getCorner(2)->z());

		y = (tets[i]->getCorner(1)->z() - tets[i]->getCorner(2)->z())*(tets[i]->getCorner(0)->x() - tets[i]->getCorner(2)->x()) -
			(tets[i]->getCorner(0)->z() - tets[i]->getCorner(2)->z())*(tets[i]->getCorner(1)->x() - tets[i]->getCorner(2)->x());

		z = (tets[i]->getCorner(1)->x() - tets[i]->getCorner(2)->x())*(tets[i]->getCorner(0)->y() - tets[i]->getCorner(2)->y()) -
			(tets[i]->getCorner(0)->x() - tets[i]->getCorner(2)->x())*(tets[i]->getCorner(1)->y() - tets[i]->getCorner(2)->y());

		x_m = (tets[i]->getCorner(0)->x() + tets[i]->getCorner(1)->x() + tets[i]->getCorner(2)->x()) / 3;
		y_m = (tets[i]->getCorner(0)->y() + tets[i]->getCorner(1)->y() + tets[i]->getCorner(2)->y()) / 3;
		z_m = (tets[i]->getCorner(0)->z() + tets[i]->getCorner(1)->z() + tets[i]->getCorner(2)->z()) / 3;

		glVertex3d(x_m, y_m, z_m);
		glVertex3d(x_m + x / 200, y_m + y / 200, z_m + z / 200);
	}
	glEnd();

}

void display()
{
	static float angle = 0.0f;

	glClear(GL_COLOR_BUFFER_BIT);
	glPushMatrix();
	glTranslatef(0, 0, -300);
	glRotatef(angle, 0.0f, 1.0f, 0.0f);

	/*
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 1000,0);
	glEnd();
	*/

	displayCorners();
	//displayVectorsB();

	glPopMatrix();
	glutSwapBuffers();

	angle += 0.02f;
}

void initGL(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

	glutInitWindowSize(1000, 1000);
	glutCreateWindow("HydroDynamics");
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glFrustum(-100, 100, -100, 100, 100, 2000);  //sets the perspective of the image. More details: opengl.gamedev.ru/doc/?func=glFrustum
	glutDisplayFunc(display);
	//glutIdleFunc(display);
	glutMainLoop();


}
#endif // DRAWMODULE_H
