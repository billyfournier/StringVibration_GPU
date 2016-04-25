//gcc graphingfxt.c -o graph -lglut -lm -lGLU -lGL
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265359 

#define L 0.67945       //Length of string in Meters
#define D 0.18415       //Laterail displacement of string in Meters
#define MASS  0.006     //Mass of string in Kilograms
#define K 38.6220      	//Spring constant of string in Newton/Meters
#define Tention 1.0 //5.3955 	//Resting tention in Newtons
#define DAMP 0.00005   	//Air resistance

#define X_WINDOW 1000
#define Y_WINDOW 700

#define X_MAX L //0.54
#define X_MIN 0.0 //0.0
#define X_SCALE 0.1

#define Y_MAX D //0.025
#define Y_MIN -D //-0.025
#define Y_SCALE 0.01

// globals
const float g_time_step = 0.00001;
const int g_number_of_segments = 1000;
const float g_time_stop = 100.0;
const int g_terms = 25;
const float g_L = L; //0.54;
const float g_d = D; //0.025;
const float g_a = 863.2596685;
const float g_damp = 4.0;

/*	Takes machine x and y which start in the upper left corner and go from zero to X_WINDOW
	left to right and form zero to Y_WINDOW top to bottom and transslates this into screen 
	points which are a -1 to 1, -1 to 1 window.
*/
float x_machine_to_x_screen(int x)
{
	return( (2.0*x)/X_WINDOW-1.0 );
}

float y_machine_to_y_screen(int y)
{
	return( -(2.0*y)/Y_WINDOW+1.0 );
}

/*	Takes machine x and y which start in the upper left corner and go from zero to X_WINDOW
	left to right and form zero to Y_WINDOW top to bottom and translates this into world 
	points which are a X_MIN to X_MAX, Y_MIN to Y_MAX window.
*/
float x_machine_to_x_world(int x)
{
	float range;
	range = X_MAX - X_MIN;
	return( (range/X_WINDOW)*x + X_MIN );
}

float y_machine_to_y_world(int y)
{
	float range;
	range = Y_MAX - Y_MIN;
	return(-((range/Y_WINDOW)*y - Y_MAX));
}

/*	Take world  points to screen points 
*/
float x_world_to_x_screen(float x)
{
	float range;
	range = X_MAX - X_MIN;
	return( -1.0 + 2.0*(x - X_MIN)/range );
}

float y_world_to_y_screen(float y)
{
	float range;
	range = Y_MAX - Y_MIN;
	return( -1.0 + 2.0*(y - Y_MIN)/range );
}

void place_axis()
{
	glColor3f(1.0,0.0,0.0);

	glBegin(GL_LINE_LOOP);
		glVertex2f(x_world_to_x_screen(X_MIN),y_world_to_y_screen(0.0));
		glVertex2f(x_world_to_x_screen(X_MAX),y_world_to_y_screen(0.0));
	glEnd();

	glBegin(GL_LINE_LOOP);
		glVertex2f(x_world_to_x_screen(0.0),y_world_to_y_screen(Y_MIN));
		glVertex2f(x_world_to_x_screen(0.0),y_world_to_y_screen(Y_MAX));
	glEnd();

	glutSwapBuffers();
}

void place_hash_marks()
{
	float x,y,dx,dy;

	glColor3f(1.0,0.0,0.0);

	dx = X_SCALE;
	dy = Y_SCALE;

	x = X_MIN;
	while(x <= X_MAX)
	{
		glBegin(GL_LINE_LOOP);
			glVertex2f(x_world_to_x_screen(x), 0.025+y_world_to_y_screen(0));
			glVertex2f(x_world_to_x_screen(x),-0.025+y_world_to_y_screen(0));
		glEnd();

		x = x + dx;
	}

	y = Y_MIN;
	while(y <= Y_MAX)
	{
		glBegin(GL_LINE_LOOP);
			glVertex2f( 0.005+x_world_to_x_screen(0),y_world_to_y_screen(y));
			glVertex2f(-0.005+x_world_to_x_screen(0),y_world_to_y_screen(y));
		glEnd();

		y = y + dy;
	}
}

float f(float x, float t)
{
	int n;
	float temp, a, Bn, nn;
	
	a = sqrt(g_a);
	
	//t = 0.0;
	
	temp = 0.0;
	for(n = 1; n <= g_terms; n++)
	{
		nn = 2.0*(float)n - 1.0;
		
		Bn = (8.0*g_d)/(nn*nn*PI*PI);
		if(n%2 == 0) Bn = -1.0*Bn;
		temp += Bn*cos((a*nn*PI*t)/g_L)*sin((nn*PI*x)/g_L);
	}
	return(temp*exp(-g_damp*t));
}

void draw_segment(float x,float t, float dx)
{
	glColor3f(1.0,0.0,0.0);
	glBegin(GL_LINE_LOOP);
		glVertex2f(x_world_to_x_screen(x),y_world_to_y_screen(f(x,t)));
		glVertex2f(x_world_to_x_screen(x+dx),y_world_to_y_screen(f(x+dx, t)));
	glEnd();
}

void graph()
{
	int i;
	float range,dx,dt;
	float t,x;
	
	range = X_MAX - X_MIN;
	dx = range/g_number_of_segments;
	dt = g_time_step;
	x = X_MIN;
	t = 0.0;
	
	printf("\ninter a character in the terminal\n");
	getchar();
	
	while(t <= g_time_stop)
	{
		glClear(GL_COLOR_BUFFER_BIT);
		//place_axis();
		//place_hash_marks();
		
		while(x <= X_MAX)
		{
			draw_segment(x,t, dx);
			x += dx;
		}
		x = X_MIN;
		t += dt;
		glutSwapBuffers();
	}
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	

	//place_axis();

	//place_hash_marks();

	graph();
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(X_WINDOW,Y_WINDOW);
	glutInitWindowPosition(0,0);
	glutCreateWindow("BOX");
	glutDisplayFunc(display);
	glutMainLoop();
}

