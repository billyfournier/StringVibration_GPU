//nvcc particleString.cu -o temp -lglut -lGL -lm
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265359 

#define X_WINDOW 1000
#define Y_WINDOW 700

#define L 0.67945       //Length of string in Meters
#define D 0.18415       //Laterail displacement of string in Meters
#define MASS  0.006     //Mass of string in Kilograms
#define K 38.6220      	//Spring constant of string in Newton/Meters
#define Tention 1.0 //5.3955 	//Resting tention in Newtons
#define N 1002          //number of bodies
#define P 1002	        //number of bodies per block
#define DAMP 0.00005   	//Air resistance

#define X_MAX (L/2.0)
#define X_MIN -(L/2.0)
#define X_SCALE 0.1

#define Y_MAX D
#define Y_MIN -D
#define Y_SCALE 0.1

#define TIME_DURATION	1000.0
#define STEP_SIZE        0.0000005
#define TIME_STEP_BETWEEN_VIEWING 100 
float *X_CPU, *Y_CPU, *VX_CPU, *VY_CPU, *AX_CPU, *AY_CPU;  //CPU pointers

float *X_GPU, *Y_GPU, *VX_GPU, *VY_GPU, *AX_GPU, *AY_GPU; //GPU pointers

dim3 dimBlock, dimGrid; //Block and Grid Dimensions


void 	SetUpCudaDevices() // Sets up the architecture for processes
	{	
		//Threads in a block
		dimBlock.x = P;
		dimBlock.y = 1;
		dimBlock.z = 1;
	
		//Blocks in a grid
		dimGrid.x = 1;
		dimGrid.y = 1;
		dimGrid.z = 1;
	}

void 	AllocateMemory()
{					
	//Allocate Device (GPU) Memory, & allocates the value of the specific pointer/array
	cudaMalloc(&X_GPU, N*sizeof(float));
	cudaMalloc(&Y_GPU, N*sizeof(float));
	cudaMalloc(&VX_GPU,N*sizeof(float));
	cudaMalloc(&VY_GPU,N*sizeof(float));
	cudaMalloc(&AX_GPU,N*sizeof(float));
	cudaMalloc(&AY_GPU,N*sizeof(float));

	//Allocate Host (CPU) Memory
	X_CPU  = (float*)malloc(N*sizeof(float)); //(float*) to prevent from being a void
	Y_CPU  = (float*)malloc(N*sizeof(float));
	VX_CPU = (float*)malloc(N*sizeof(float));
	VY_CPU = (float*)malloc(N*sizeof(float));
	AX_CPU = (float*)malloc(N*sizeof(float));
	AY_CPU = (float*)malloc(N*sizeof(float));
}

float x_machine_to_x_screen(int x)
{
	return( (2.0*x)/X_WINDOW-1.0 );
}

float y_machine_to_y_screen(int y)
{
	return( -(2.0*y)/Y_WINDOW+1.0 );
}

/*	Takes machine x and y which start in the upper left corner and go from zero to X_WINDOW
	left to right and form zero to Y_WINDOW top to bottom and transslates this into world 
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

void 	draw_spring(float *x, float *y)
	{
		int i;

		glPointSize(1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		glColor3f(1.0,1.0,0.0);
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(x[0]),y_world_to_y_screen(y[0]));
		glEnd();

		glColor3f(1.0,1.0,1.0);
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(x[N-1]),y_world_to_y_screen(y[N-1]));
		glEnd();
			
		glColor3f(1.0,0.0,0.0);
		for(i = 1; i < N-1; i++)
		{			
			glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(x[i]),y_world_to_y_screen(y[i]));
			glEnd();
			
		}
		glFlush();
	}

__global__ void Findforce(float *x, float *y, float *vx, float *vy, float *ax, float *ay, float l, float mass, float k) // This is the kernel, it is the function that is being fed to the GPU.
{	
	int id = threadIdx.x;
	float dx, dy, d2, d, f;

	if(0 < threadIdx.x && threadIdx.x < P-1)
	{
		dx = x[id-1]-x[id];
		dy = y[id-1]-y[id];
		d2 = dx*dx + dy*dy;
		d  = sqrt(d2);
		f = k*(d-l) + Tention;

		ax[id] += (f*dx/d)/mass;
		ay[id] += (f*dy/d)/mass;
			
			
		dx = x[id+1]-x[id];
		dy = y[id+1]-y[id];
		d2 = dx*dx + dy*dy;
		d  = sqrt(d2);
		f = k*(d-l) + Tention;		
		ax[id] += (f*dx/d)/mass;
		ay[id] += (f*dy/d)/mass;
		
		ax[id] += (-DAMP*vx[id])/mass;
		ay[id] += (-DAMP*vy[id])/mass;
	}
	__syncthreads();
}

void n_body()
{
	float l, time, dt;
	int draw_count,i;
	float particleMass, particleK;

	SetUpCudaDevices();
	AllocateMemory();

	l = L/(N-1);

	time = 0.0;
	dt = STEP_SIZE;
	draw_count = 0;
	
	particleMass = (MASS/L)/(float)N;
	particleK = K*(float(N-1));
	//particleK = K;

	X_CPU[0] = -L/2.0;
	Y_CPU[0] = 0.0;
	X_CPU[N-1] = L/2.0;
	Y_CPU[N-1] = 0.0;
	
	for(i=1; i<(N-1); i++)
	{
		X_CPU[i] = -L/2.0 + l*(i);
		if(X_CPU[i] <= 0.0) Y_CPU[i] = (2.0*D/L)*(X_CPU[i]+L/2.0);
		if(X_CPU[i] >  0.0) Y_CPU[i] = D - (2.0*D/L)*(X_CPU[i]);
		VX_CPU[i] = 0.0;
		VY_CPU[i] = 0.0;
	}
	draw_spring(X_CPU, Y_CPU);
	
	printf("\ninter a character in the terminal\n");
	getchar();
	
	while(time < TIME_DURATION)
	{
		for(i=0; i<N; i++) 
		{
			AX_CPU[i] = 0.0;
			AY_CPU[i] = 0.0;
		}	
			
		//Copy Memory from CPU to GPU		
		cudaMemcpyAsync(X_GPU,   X_CPU, N*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(Y_GPU,   Y_CPU, N*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(VX_GPU, VX_CPU, N*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(VY_GPU, VY_CPU, N*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(AX_GPU, AX_CPU, N*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(AY_GPU, AY_CPU, N*sizeof(float), cudaMemcpyHostToDevice);
	
		//Launch Kernel	
		Findforce<<<dimGrid, dimBlock>>>(X_GPU, Y_GPU, VX_GPU, VY_GPU, AX_GPU, AY_GPU, l, particleMass, particleK);
	
		//Copy Memory from GPU to CPU	
		cudaMemcpyAsync(AY_CPU, AY_GPU, N*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpyAsync(AX_GPU, AX_CPU, N*sizeof(float), cudaMemcpyDeviceToHost);
	
		for(i=1; i < (N-1); i++)
		{
			VX_CPU[i] += AX_CPU[i]*dt;
			VY_CPU[i] += AY_CPU[i]*dt;
			X_CPU[i]  += VX_CPU[i]*dt;
			Y_CPU[i]  += VY_CPU[i]*dt;
		}

		if(draw_count == TIME_STEP_BETWEEN_VIEWING)
		{
			draw_spring(X_CPU, Y_CPU);
			draw_count = 0;
		}

		time = time + dt;
		draw_count++;
	}
}

void display()
{
	//glClear(GL_COLOR_BUFFER_BIT);
	//glFlush();
	
	n_body();
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitWindowSize(X_WINDOW,Y_WINDOW);
	glutInitWindowPosition(0,0);
	glutCreateWindow("BOX");
	glutDisplayFunc(display);
	glutMainLoop();
}




