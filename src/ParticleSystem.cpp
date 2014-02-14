#include "ParticleSystem.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "vec2.h"

#define PI 3.1415926535f

#define MAX_PARTICLES 16384
#define MAX_BORDER_PARTICLES 5000
#define averageParticles 20

#define BORDER_LEFT -0.8
#define BORDER_RIGHT 0.8
#define BORDER_TOP 0.8
#define BORDER_BOTTOM -0.8
#define SIM_WIDTH 1.6
#define SIM_HEIGHT 1.6

#define SIM_SCALE (SIM_WIDTH/2)
#define TIME_STEP 0.0005 //s
#define INTERACTION_RADIUS 0.005 //m

#define INTERACTION_RADIUS2 (INTERACTION_RADIUS*INTERACTION_RADIUS)

#define CELL_SIZE INTERACTION_RADIUS
#define GRID_WIDTH (SIM_WIDTH/CELL_SIZE)
#define GRID_HEIGHT (SIM_HEIGHT/CELL_SIZE)

#define EPSILON	0.0000001f			//for collision detection


#define STIFFNESS 12.0 //
#define VISCOSITY 3.5 // pascal-seconds
//#define PARTICLE_MASS 0.00020543 //kg


#define REST_DENSITY 1000.0 //kg / m^3
#define VEL_LIMIT 200.0 //velocity limit (m/s)
#define PARTICLE_RADIUS 0.004 // m
#define EXT_DAMP 512.0

float PARTICLE_MASS = 0.0f;// 0.000159963f;

const float W_DEFAULT = 315.0f / (64.0f * 3.141592 * pow( INTERACTION_RADIUS, 9) );
const float	W_GRAD_PRESSURE = -45.0f / (3.141592 * pow( INTERACTION_RADIUS, 6) );
const float W_LAPLACIAN_VISCOSITY = 45.0f / (3.141592 * pow( INTERACTION_RADIUS, 6) );

ParticleSystem::ParticleSystem(void)
{
	particles.resize(MAX_PARTICLES + MAX_BORDER_PARTICLES);
	border_particles.resize(MAX_BORDER_PARTICLES);
	grid.resize(GRID_WIDTH*GRID_HEIGHT);
	borderParticleCount = 0;
	createParticleField();
	//createBorderParticles();
	updateGrid();
	calculateMass();
	std::cout << "Particle mass: " << PARTICLE_MASS << std::endl;
	advance_call = 0;
	draw_counter = 0;
	particleCount = 0;
	emitStep = 0;
	borderStep = 0;
}

void ParticleSystem::reloadParticleSystem()
{
	memset(&particles[0], 0, particles.size()*sizeof(particles[0]));
	createParticleField();
	updateGrid();
	advance_call = 0;
	draw_counter = 0;
}


void ParticleSystem::reloadParticleSystem(char c)
{
	memset(&particles[0], 0, particles.size()*sizeof(particles[0]));
	particles.resize(MAX_PARTICLES + borderParticleCount);
	int rowcolSize = sqrt(MAX_PARTICLES);
	keyPressed = c;
	switch (c)
	{
	case '1':
		particleCount = MAX_PARTICLES;
		borderParticleCount = 0;
		for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
		{
			float stepLength = 0.0038; //(BORDER_RIGHT - BORDER_LEFT)/(3*rowcolSize);
			for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
			{
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_LEFT+0.002 + particleIndexCol*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_BOTTOM+0.002 + particleIndexRow*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
			}
		}
		break;

	case '2':
		particleCount = MAX_PARTICLES;
		borderParticleCount = 0;
		for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
		{
			float stepLength = 0.0038; //(BORDER_RIGHT - BORDER_LEFT)/(3*rowcolSize);
			for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
			{
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_LEFT+0.002 + particleIndexCol*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_BOTTOM+0.2 + particleIndexRow*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
			}
		}
		break;

	case '3':
		particleCount = MAX_PARTICLES;
		borderParticleCount = 0;
		for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
		{
			float stepLength = 0.0038; //(BORDER_RIGHT - BORDER_LEFT)/(3*rowcolSize);
			for(int particleIndexCol = 0; particleIndexCol < (rowcolSize/2); ++particleIndexCol)
			{
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_LEFT+0.002 + particleIndexCol*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_BOTTOM+0.2 + particleIndexRow*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
			}
		}
		for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
		{
			float stepLength = 0.0038; //(BORDER_RIGHT - BORDER_LEFT)/(3*rowcolSize);
			for(int particleIndexCol = (rowcolSize/2); particleIndexCol < rowcolSize; ++particleIndexCol)
			{
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_RIGHT-stepLength*(rowcolSize) + particleIndexCol*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_BOTTOM+0.2 + particleIndexRow*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
			}
		}
		break;

	case '4':
		particleCount = MAX_PARTICLES;
		borderParticleCount = 0;
		for(int particleIndexRow = 0; particleIndexRow < (rowcolSize/2); ++particleIndexRow)
		{
			float stepLength = 0.0038; //(BORDER_RIGHT - BORDER_LEFT)/(3*rowcolSize);
			for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
			{
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_LEFT+0.01 + particleIndexCol*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_TOP-stepLength*rowcolSize + particleIndexRow*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].velocity.x = 1;
				particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
			}
		}
		for(int particleIndexRow = (rowcolSize/2); particleIndexRow < rowcolSize; ++particleIndexRow)
		{
			float stepLength = 0.0038; //(BORDER_RIGHT - BORDER_LEFT)/(3*rowcolSize);
			for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
			{
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_RIGHT-stepLength*(rowcolSize) - 0.005 + particleIndexCol*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_TOP-stepLength*rowcolSize - 0.02 + particleIndexRow*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].velocity.x = -1;
				particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
			}
		}
		break;
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
		particleCount=0;
		emitStep = 0;
		particles.resize(0);
		break;
	case '0':
		particleCount = MAX_PARTICLES;
		borderParticleCount = 0;
		for(int i = 0; i < MAX_PARTICLES; i++)
		{
			particles[i].position.x =Random(-0.02,0.02);
			particles[i].position.y = Random(-0.02,0.02);
			particles[i].is_static = false;
		}
		break;
	default:
		particleCount = MAX_PARTICLES;
		for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
		{
			float stepLength = 0.0038; //(BORDER_RIGHT - BORDER_LEFT)/(3*rowcolSize);
			for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
			{
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_LEFT+0.002 + particleIndexCol*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_BOTTOM+0.002 + particleIndexRow*stepLength;
				particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
			}
		}
		break;
	}
	updateGrid();
	advance_call = 0;
	draw_counter = 0;
}

std::vector<vec2> ParticleSystem::getParticleCoordinates()
{
	std::vector<vec2> coordinateVector(particles.size());
	for(int i = 0; i < particles.size(); ++i)
	{
		coordinateVector[i] = particles[i].position/SIM_SCALE;
	}

	/*for(int i=0; i < border_particles.size(); ++i)
	{
		coordinateVector[particles.size()+i] = border_particles[i].position/SIM_SCALE;
	}*/

	return coordinateVector;
}

std::vector<vec2> ParticleSystem::getParticleCoordinatesBorder()
{
	std::vector<vec2> coordinateVector(border_particles.size());
	for(int i = 0; i < border_particles.size(); ++i)
	{
		coordinateVector[i] = border_particles[i].position/SIM_SCALE;
	}

	return coordinateVector;
}

std::vector<vec2> ParticleSystem::getParticleCoordinatesNeighbours()
{
	int particle_index = (advance_call/20)%MAX_PARTICLES;
	std::vector<vec2> coordinateVector(particles[particle_index].neighbour_count);
	for(int i = 0; i < particles[particle_index].neighbour_count; ++i)
		coordinateVector[i] = particles[particle_index].neighbours[i]->position/SIM_SCALE;

	return coordinateVector;
}

std::vector<vec2> ParticleSystem::getParticleCoordinatesGrid()
{
	int grid_index = (advance_call)%grid.size();
	std::vector<vec2> coordinateVector;
	int j=0;
	for (particle* ppj=grid[grid_index]; ppj != 0; ppj=ppj->next)
		coordinateVector.push_back(ppj->position/SIM_SCALE);

	return coordinateVector;
}

std::vector<vec2> ParticleSystem::getParticleCoordinatesPressure(float dir, float limit)
{
	std::vector<vec2> coordinateVector;
	for(particle& pi : particles)
		if(pi.pressure*dir >= limit*dir)
			coordinateVector.push_back(pi.position/SIM_SCALE);

	return coordinateVector;
}


#define OFFSET BORDER_LEFT + 0.1
void ParticleSystem::createParticleField()
{
	particles.resize(MAX_PARTICLES+borderParticleCount);

	int rowcolSize = sqrt(MAX_PARTICLES);
	particleCount = MAX_PARTICLES;
	for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
	{
		float stepLength = 0.0038; //(BORDER_RIGHT - BORDER_LEFT)/(3*rowcolSize);
		for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
		{
			particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_LEFT+0.002 + particleIndexCol*stepLength;
			particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_BOTTOM+0.002 + particleIndexRow*stepLength;
			particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
			particles[particleIndexCol + rowcolSize*particleIndexRow].mass = PARTICLE_MASS;
		}
	}
}

void ParticleSystem::createBorderParticles()
{
	border_particles.resize(MAX_BORDER_PARTICLES);
	float stepLength = (BORDER_RIGHT - BORDER_LEFT)/(MAX_BORDER_PARTICLES*SIM_SCALE);
	for( int i = 0; i < MAX_BORDER_PARTICLES; ++i)
	{
		border_particles[i].pressure = 1000;
		border_particles[i].density = 1000;
		border_particles[i].position.x = BORDER_RIGHT;//BORDER_LEFT + i*stepLength;
		border_particles[i].position.y = BORDER_TOP;//BORDER_BOTTOM+0.1;
	}
}


void ParticleSystem::updateGrid()
{
	memset(&grid[0], 0, grid.size()*sizeof(grid[0]));
	for(particle& pi : particles)
	{
		int x = (-(BORDER_LEFT) + pi.position.x)/CELL_SIZE;
		int y = (-(BORDER_BOTTOM) + pi.position.y)/CELL_SIZE;

		if (x < 1)
			x = 1;
		else if (x > GRID_WIDTH-2)
			x = GRID_WIDTH-2;

		if (y < 1)
			y = 1;
		else if (y > GRID_HEIGHT-2)
			y = GRID_HEIGHT-2;

		pi.next = grid[x+y*GRID_WIDTH];
		grid[x+y*GRID_WIDTH] = &pi;

		pi.grid_x = x;
		pi.grid_y = y;
	}
}

void ParticleSystem::updateNeighbours()
{
	#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < particles.size(); i++)
	{
		particle& pi = particles[i];
		pi.neighbour_count = 0;
		size_t gi = pi.grid_x;
		size_t gj = pi.grid_y*GRID_WIDTH;

		

		/*for(int j=0; j<border_particles.size(); j++)
		{
			particle& pj = border_particles[j];
			vec2 distance_vector = (pi.position - pj.position);
			float distance2 = dot(distance_vector,distance_vector);
			if(distance2 < INTERACTION_RADIUS2)
			{
				if(pi.neighbour_count < MAX_NEIGHBOURS)
				{
					pi.neighbours[pi.neighbour_count] = &pj;
					pi.neighbour_distance[pi.neighbour_count] = sqrt(distance2);
					++pi.neighbour_count;
					//std::cout << "I'm on the border" << std::endl;
				}
			}
		}*/
		
		//loop over adjacent cells
		#pragma omp parallel for
		for (int ni=gi-1; ni<=gi+1; ++ni)
		{
			for (int nj=gj-GRID_WIDTH; nj<=gj+GRID_WIDTH; nj+=GRID_WIDTH)
			{
				//loop over neighbours
				for(particle* ppj=grid[ni+nj]; ppj != 0; ppj=ppj->next)
				{
					if(&pi != ppj)
					{
						vec2 distance_vector = (pi.position - ppj->position);
						float distance2 = dot(distance_vector,distance_vector);
						if(distance2 < INTERACTION_RADIUS2)
						{
							if(pi.neighbour_count < MAX_NEIGHBOURS)
							{
								pi.neighbours[pi.neighbour_count] = ppj;
								pi.neighbour_distance[pi.neighbour_count] = sqrt(distance2);
								++pi.neighbour_count;
							}
						}
					}
				}
			}
		}
	}
}

void ParticleSystem::calculateMass()
{
	float density = 0.0f;

	for(particle& pi : particles)
	{
		size_t gi = pi.grid_x;
		size_t gj = pi.grid_y*GRID_WIDTH;
		//loop over adjacent cells
		for (int ni=gi-1; ni<=gi+1; ++ni)
		{
			for (int nj=gj-GRID_WIDTH; nj<=gj+GRID_WIDTH; nj+=GRID_WIDTH)
			{
				//loop over neighbours
				for (particle* ppj=grid[ni+nj]; ppj != 0; ppj=ppj->next)
				{
					if(&pi != ppj)
					{
						vec2 distance_vector = (pi.position - ppj->position);
						float distance2 = dot(distance_vector,distance_vector);
						if(distance2 < INTERACTION_RADIUS2)
						{
								density += W_DEFAULT * (INTERACTION_RADIUS2 - distance2)*(INTERACTION_RADIUS2 - distance2)*(INTERACTION_RADIUS2 - distance2);
						}
					}
				}
			}
		}
	}

	float dA = density/MAX_PARTICLES;
	PARTICLE_MASS = (dA*REST_DENSITY)/(dA*dA);


}

void ParticleSystem::calculatePressure()
{
	//float c;
	#pragma omp parallel for schedule(guided)
	for(int j = 0; j < MAX_PARTICLES; j++)
	{
		particle& pi = particles[j];
		float sum = 0.0;
		int i;	

		//#pragma omp parallel for schedule(static)
		for(i = 0; i < pi.neighbour_count; i++)
		{
			particle& pj = *(pi.neighbours[i]);
			vec2 distance_vector = (pi.position - pj.position);
			float distance2 = dot(distance_vector,distance_vector);
			if(distance2 < INTERACTION_RADIUS2)
			{
				float mass;
				if(pj.is_static)
				{
					mass = pj.mass;
				}
				else
				{
					mass = PARTICLE_MASS;
				}

				float c = (INTERACTION_RADIUS2 - distance2);
				//#pragma omp atomic
				sum += c*c*c*mass;
			}
		}
		pi.density = sum * W_DEFAULT;
		pi.pressure = (pi.density - REST_DENSITY) * STIFFNESS;
		pi.density = 1.0f/pi.density;
	}
}

void ParticleSystem::calculateSPHForce()
{
	#pragma omp parallel for schedule(guided)
	for(int j = 0; j < MAX_PARTICLES; j++)
	{
		particle& pi = particles[j];
		if(pi.is_static) continue;
		vec2 force;
		float forcex = 0;
		float forcey = 0;
		int i;
		//#pragma omp parallel for 
		for(i=0; i<pi.neighbour_count; i++)
		{
			particle& pj = *(pi.neighbours[i]);
			vec2 distance_vector = (pi.position - pj.position);
			float distance2 = dot(distance_vector,distance_vector);
			float c = ( INTERACTION_RADIUS - pi.neighbour_distance[i] );
			float pressure_term = -0.5f * c * W_GRAD_PRESSURE * (pi.pressure + pj.pressure) / pi.neighbour_distance[i];
			float density_term = c * pi.density * pj.density;
			float viscosity_term = W_LAPLACIAN_VISCOSITY * VISCOSITY;
			//#pragma omp atomic
			forcex += ( pressure_term * distance_vector.x + viscosity_term * (pj.velocity_eval.x - pi.velocity_eval.x) ) * density_term;
			//#pragma omp atomic
			forcey += ( pressure_term * distance_vector.y + viscosity_term * (pj.velocity_eval.y - pi.velocity_eval.y) ) * density_term;
			if(pj.is_static)
			{
				forcex = 0;//pressure_term;
				forcey = 0;//pressure_term;
				float A = atan2(pi.position.y - pj.position.y,pi.position.x-pj.position.x);
				float vx = cos(pi.velocity.x-A);
				float vy = cos(pi.velocity.y-A);

				float fx = vx*(pi.mass - pj.mass);
				float fy = vy*(pi.mass - pj.mass);
				pi.velocity.x = 0.2*(fx+A);
				pi.velocity.y = 0.2*(fy+A);


			}

			//force += ( pressure_term * distance_vector + viscosity_term * (pj.velocity_eval - pi.velocity_eval) ) * density_term;
		}
		pi.force.x = forcex;
		pi.force.y = forcey;

	}
}

void ParticleSystem::moveParticles()
{

	float speed,diff,adj;

	float current, u, v, cp, n,d;

	float SL2 = VEL_LIMIT*VEL_LIMIT;
	vec2 norm;

	#pragma omp parallel for schedule(guided)
	for(int i = 0; i < particles.size(); i++)
	{
		particle& pi = particles[i];
		if(pi.is_static) continue;
		vec2 acceleration = pi.force*PARTICLE_MASS;

		speed = acceleration.x*acceleration.x + acceleration.y*acceleration.y;
		if(speed > SL2)
		{
			acceleration *= VEL_LIMIT / sqrt(speed);
		}

		acceleration.y += -9.81;
		
		diff = 2 * PARTICLE_RADIUS - (pi.position.y - BORDER_LEFT);
		if( diff > EPSILON)
		{
			norm.x = 0.0f;
			norm.y = 1.0f;

			adj = STIFFNESS * diff - EXT_DAMP * dot(norm, pi.velocity_eval);
			acceleration.x += adj * norm.x;
			acceleration.y += adj * norm.y;
		}

		diff = 2 * PARTICLE_RADIUS - (BORDER_RIGHT - pi.position.y);
		if( diff > EPSILON)
		{
			norm.x = 0.0f;
			norm.y = -1.0f;

			adj = STIFFNESS * diff - EXT_DAMP * dot(norm, pi.velocity_eval);
			acceleration.x += adj * norm.x;
			acceleration.y += adj * norm.y;
		}

		diff = 2 * PARTICLE_RADIUS - (pi.position.x - BORDER_LEFT);
		if( diff > EPSILON)
		{
			norm.x = 1.0f;
			norm.y = 0.0f;

			adj = STIFFNESS * diff - EXT_DAMP * dot(norm, pi.velocity_eval);
			acceleration.x += adj * norm.x;
			acceleration.y += adj * norm.y;
		}

		diff = 2 * PARTICLE_RADIUS - (BORDER_RIGHT - pi.position.x);
		if( diff > EPSILON)
		{
			norm.x = -1.0f;
			norm.y = 0.0f;

			adj = STIFFNESS * diff - EXT_DAMP * dot(norm, pi.velocity_eval);
			acceleration.x += adj * norm.x;
			acceleration.y += adj * norm.y;
		}

		if(pi.position.y < BORDER_BOTTOM)
			pi.position.y = BORDER_BOTTOM;

		if(pi.position.x < BORDER_LEFT)
			pi.position.x = BORDER_LEFT;

		if(pi.position.x > BORDER_RIGHT)
			pi.position.x = BORDER_RIGHT;
		

		//leapfrog integration
		vec2 velocity_next = pi.velocity + acceleration*TIME_STEP;
		pi.velocity_eval = (2.0*pi.velocity + acceleration*TIME_STEP)*0.5;
		pi.velocity = velocity_next;
		velocity_next *= TIME_STEP;
		pi.position += velocity_next;


		/*if(pi.position.y < BORDER_BOTTOM)
		{
			current = pi.position.y;
			u = pi.velocity.x;
			v = pi.velocity.y;
			cp = BORDER_BOTTOM;

			d = sqrt((cp-current)*(cp-current));
			n = 1;
			pi.position.y = cp + d*n;
			pi.velocity_eval.y = 0.0f; //v - (1 + EXT_DAMP*(d/(TIME_STEP*sqrt(u*u+v*v))))*(v*n)*n;
			pi.velocity_eval.x = u;
		}

		if(pi.position.x < BORDER_LEFT)
		{
			current = pi.position.x;
			u = pi.velocity.x;
			v = pi.velocity.y;
			cp = BORDER_LEFT;

			d = sqrt((cp-current)*(cp-current));
			n = 1;
			pi.position.x = cp + d*n;
			pi.velocity_eval.x = 0.0f; //u - (1 + EXT_DAMP*(d/(TIME_STEP*sqrt(u*u+v*v))))*(u*n)*n;
			pi.velocity_eval.y = v;
		}

		if(pi.position.x > BORDER_RIGHT)
		{
			current = pi.position.x;
			u = pi.velocity.x;
			v = pi.velocity.y;
			cp = BORDER_RIGHT;

			d = sqrt((cp-current)*(cp-current));
			n = 1;
			pi.position.x = cp - d*n;
			pi.velocity_eval.x = 0.0f; //u - (1 + EXT_DAMP*(d/(TIME_STEP*sqrt(u*u+v*v))))*(u*n)*n;
			pi.velocity_eval.y = v;
		}*/

	}

}

void ParticleSystem::advance()
{
	if(keyPressed == '5' || keyPressed == '6' || keyPressed == '7'
		|| keyPressed == '8')
	{
		if(emitStep == 0)
		{
			EmitParticles();
			emitStep=6;
		}
		emitStep--;
	}
	if(rButtonPressed)
	{
		if(borderStep == 0)
		{
			setBorderParticles();
			borderStep = 10;
		}
		borderStep--;
	}
	updateGrid();
	updateNeighbours();
	calculatePressure();
	calculateSPHForce();
	moveParticles();
	advance_call++;
}

void ParticleSystem::drawParticle(float x, float y, bool is_static)
{
	particle p;
	p.position = vec2(x,y);
	p.is_static = is_static;

	if(particles.size() < MAX_PARTICLES + borderParticleCount)
	{
		particles.push_back(p);
	}
	else
	{
		particles[draw_counter%particles.size()] = p;
	}
}

void ParticleSystem::EmitParticles()
{
	int particleAddCount = 10;
	float stepLength = 0.0038;
	float totParticleCount = particleCount + borderParticleCount;
	if(keyPressed == '5')
	{
		if(particleCount+particleAddCount < MAX_PARTICLES){
			particles.resize(totParticleCount+particleAddCount);
			for(int i=0; i<particleAddCount; i++)
			{
				particles[totParticleCount+i].position.x = BORDER_LEFT+0.01+i*(stepLength/2);
				particles[totParticleCount+i].position.y = 0.5-i*stepLength;
				particles[totParticleCount+i].velocity.x = 50*particleAddCount*stepLength*Random(0.9,1.1);
				particles[totParticleCount+i].velocity.y = 50*particleAddCount*(stepLength/2)*Random(0.9,1.1);
			}
			totParticleCount += particleAddCount;
		}
	}
	if(keyPressed == '6')
	{
		if(particleCount+2*particleAddCount < MAX_PARTICLES)
		{
			particles.resize(totParticleCount+2*particleAddCount);
			for(int i=0; i<particleAddCount; i++)
			{
				particles[totParticleCount+i].position.x = BORDER_LEFT+0.01+i*(stepLength/2);
				particles[totParticleCount+i].position.y = 0.5-i*stepLength;
				particles[totParticleCount+i].velocity.x = 50*particleAddCount*stepLength*Random(0.9,1.1);
				particles[totParticleCount+i].velocity.y = 50*particleAddCount*(stepLength/2)*Random(0.9,1.1);
				
				particles[totParticleCount+i+particleAddCount].position.x = BORDER_RIGHT-0.01-i*(stepLength/2);
				particles[totParticleCount+i+particleAddCount].position.y = 0.5-i*stepLength;
				particles[totParticleCount+i+particleAddCount].velocity.x = -50*particleAddCount*stepLength*Random(0.9,1.1);
				particles[totParticleCount+i+particleAddCount].velocity.y = 50*particleAddCount*(stepLength/2)*Random(0.9,1.1);
				
			}
			particleCount += 2*particleAddCount;
		}
	}

	if(keyPressed == '7')
	{
		if(particleCount+2*particleAddCount < MAX_PARTICLES)
		{
			particles.resize(totParticleCount+2*particleAddCount);
			for(int i=0; i<particleAddCount; i++)
			{
				particles[totParticleCount+i].position.x = BORDER_LEFT+0.01+i*(stepLength/2);
				particles[totParticleCount+i].position.y = -0.1-i*stepLength;
				particles[totParticleCount+i].velocity.x = 60*particleAddCount*stepLength*Random(0.9,1.1);
				particles[totParticleCount+i].velocity.y = 60*particleAddCount*(stepLength/2)*Random(0.9,1.1);
			
				particles[totParticleCount+i+particleAddCount].position.x = BORDER_LEFT+0.01+i*(stepLength/2);
				particles[totParticleCount+i+particleAddCount].position.y = 0.1-i*stepLength;
				particles[totParticleCount+i+particleAddCount].velocity.x = 50*particleAddCount*stepLength*Random(0.9,1.1);
				particles[totParticleCount+i+particleAddCount].velocity.y = 50*particleAddCount*(stepLength/2)*Random(0.9,1.1);
			
			}
			particleCount += 2*particleAddCount;
		}
	}
	if(keyPressed == '8')
		if(lButtonPressed)
		{
			if(particleCount+particleAddCount < MAX_PARTICLES){
			particles.resize(totParticleCount+particleAddCount);
			for(int i=0; i<particleAddCount; i++)
			{
				particles[totParticleCount+i].position.x = mouseX+Random(-0.04,0.04);
				particles[totParticleCount+i].position.y = mouseY+Random(-0.04,0.04);
				/*particles[particleCount+i].velocity.x = 50*particleAddCount*stepLength*Random(0.9,1.1);*/
				particles[totParticleCount+i].velocity.y = -0.5;//50*particleAddCount*(stepLength/2)*Random(0.9,1.1);
				
			
			}
			particleCount += particleAddCount;
		}
		//std::cout << "Mouse is at " << mouseX << " and " << mouseY << std::endl;
		}
}

void ParticleSystem::updateMouseState(float x, float y, bool lpressed, bool rpressed)
{
		mouseX = x*(BORDER_RIGHT-BORDER_LEFT)-BORDER_RIGHT;
		mouseY = BORDER_TOP-y*(BORDER_TOP - BORDER_BOTTOM);
		lButtonPressed = lpressed;
		rButtonPressed=rpressed;
}

void ParticleSystem::setBorderParticles(){
	if(borderParticleCount+10 < MAX_BORDER_PARTICLES)
	{
		float totParticleCount = particleCount + borderParticleCount;
		particles.resize(totParticleCount+10);
		float stepLength = 0.0038;
		particles[totParticleCount].position.x = mouseX;
		particles[totParticleCount].position.y = mouseY;

		particles[totParticleCount+1].position.x = mouseX;
		particles[totParticleCount+1].position.y = mouseY;

		particles[totParticleCount+2].position.x = mouseX+stepLength;
		particles[totParticleCount+2].position.y = mouseY;

		particles[totParticleCount+3].position.x = mouseX+stepLength;
		particles[totParticleCount+3].position.y = mouseY+stepLength;

		particles[totParticleCount+4].position.x = mouseX+stepLength;
		particles[totParticleCount+4].position.y = mouseY-stepLength;

		particles[totParticleCount+5].position.x = mouseX-stepLength;
		particles[totParticleCount+5].position.y = mouseY;

		particles[totParticleCount+6].position.x = mouseX-stepLength;
		particles[totParticleCount+6].position.y = mouseY+stepLength;

		particles[totParticleCount+7].position.x = mouseX-stepLength;
		particles[totParticleCount+7].position.y = mouseY-stepLength;

		particles[totParticleCount+8].position.x = mouseX;
		particles[totParticleCount+8].position.y = mouseY+stepLength;

		particles[totParticleCount+9].position.x = mouseX;
		particles[totParticleCount+9].position.y = mouseY-stepLength;

		for(int i=0; i<10; i++){
			particles[totParticleCount+i].density = 1000;
			particles[totParticleCount+i].pressure = STIFFNESS*1000;
			particles[totParticleCount+i].is_static = true;
			particles[totParticleCount+i].mass = PARTICLE_MASS;
			particles[totParticleCount+i].velocity_eval.x = 0.0f;
			particles[totParticleCount+i].velocity_eval.y = 0.0f;
			borderParticleCount++;
		}
	}
	
}