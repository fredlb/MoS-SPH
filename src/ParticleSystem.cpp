#include "ParticleSystem.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "vec2.h"

#define BORDER_LEFT -1.0
#define BORDER_RIGHT 1.0
#define BORDER_TOP 1.0
#define BORDER_BOTTOM -1.0
#define SIM_WIDTH 2.0
#define SIM_HEIGHT 2.0

#define SIM_SCALE 0.05
#define TIME_STEP 0.0003 //s
#define INTERACTION_RADIUS 0.01 //m
#define INTERACTION_RADIUS2 INTERACTION_RADIUS*INTERACTION_RADIUS

#define CELL_SIZE (INTERACTION_RADIUS/SIM_SCALE)
#define GRID_WIDTH (SIM_WIDTH/CELL_SIZE)
#define GRID_HEIGHT (SIM_HEIGHT/CELL_SIZE)

#define PI 3.1415926535f
#define MAX_PARTICLES 1024

#define STIFFNESS 0.50 //
#define VISCOSITY 0.2 // pascal-seconds
#define PARTICLE_MASS 0.00020543 //kg
#define REST_DENSITY 600.0 //kg / m^3

const float W_DEFAULT = 315.0f / (64.0f * 3.141592 * pow( INTERACTION_RADIUS, 9) );
const float	W_GRAD_PRESSURE = -45.0f / (3.141592 * pow( INTERACTION_RADIUS, 6) );
const float W_LAPLACIAN_VISCOSITY = 45.0f / (3.141592 * pow( INTERACTION_RADIUS, 6) );

ParticleSystem::ParticleSystem(void)
{
	particles.reserve(MAX_PARTICLES);
	grid.resize(GRID_WIDTH*GRID_HEIGHT);
	createParticleField();
	advance_call = 0;
}

std::vector<vec2> ParticleSystem::getParticleCoordinates()
{
	
	//standard
	std::vector<vec2> coordinateVector(particles.size());
	for(int i = 0; i < particles.size(); ++i)
	{
		coordinateVector[i] = particles[i].position;
	}
	
	
	//grid test
	/*
	int grid_index = (advance_call/10)%grid.size();
	std::vector<vec2> coordinateVector;
	int j=0;
	for (particle* ppj=grid[grid_index]; ppj != 0; ppj=ppj->next)
	{
		coordinateVector.push_back(ppj->position);
	}
	*/
	/*
	//neighbour test
	int particle_index = (advance_call/10)%MAX_PARTICLES;
	std::vector<vec2> coordinateVector(particles[particle_index].neighbour_count);
	for(int i = 0; i < particles[particle_index].neighbour_count; ++i)
	{
		coordinateVector[i] = particles[particle_index].neighbours[i]->position;
	}
	*/
	return coordinateVector;
	

}

#define OFFSET BORDER_LEFT + 0.1
void ParticleSystem::createParticleField()
{
	particles.resize(MAX_PARTICLES);

	int rowcolSize = sqrt(MAX_PARTICLES);
	for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
	{
		float stepLength = 1.5f/rowcolSize;
		for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
		{
			particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = -0.9 + particleIndexCol*stepLength;
			particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = -0.9 + particleIndexRow*stepLength;
		}
	}
}

void ParticleSystem::updateGrid()
{
	memset(&grid[0], 0, grid.size()*sizeof(grid[0]));
	for(particle& pi : particles)
	{
		int x = (1 + pi.position.x)/CELL_SIZE;
		int y = (1 + pi.position.y)/CELL_SIZE;

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
	for(particle& pi : particles)
	{
		pi.neighbour_count = 0;
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
						vec2 distance_vector = (pi.position - ppj->position)*SIM_SCALE;
						float distance2 = dot(distance_vector,distance_vector);
						if(distance2 < INTERACTION_RADIUS2)
						{
							if(pi.neighbour_count < MAX_NEIGHBOURS)
							{
								pi.neighbours[pi.neighbour_count] = ppj;
								pi.neighbour_distance[pi.neighbour_count] = sqrt(distance2);
								++pi.neighbour_count;
/*näsa*/					}
	/***/				}
	  /**/			}
	/**/		}
  /**/		}
/**/	}
	}
}

void ParticleSystem::calculatePressure()
{
	float sum, c;
	for(particle& pi : particles)
	{
		sum = 0.0;
		for(int i=0; i < pi.neighbour_count; i++)
		{
			particle& pj = *(pi.neighbours[i]);
			vec2 distance_vector = (pi.position - pj.position)*SIM_SCALE;
			float distance2 = dot(distance_vector,distance_vector);
			if(distance2 < INTERACTION_RADIUS2)
			{
				c = (INTERACTION_RADIUS2 - distance2);
				sum += c*c*c;
			}
		}
		pi.density = sum * PARTICLE_MASS * W_DEFAULT;
		pi.pressure = (pi.density - REST_DENSITY) * STIFFNESS;
		pi.density = 1.0f/pi.density;
	}
}

void ParticleSystem::calculateSPHForce()
{
	for(particle& pi : particles)
	{
		vec2 force;
		for(int i=0; i<pi.neighbour_count; i++)
		{
			particle& pj = *(pi.neighbours[i]);
			vec2 distance_vector = (pi.position - pj.position)*SIM_SCALE;
			float distance2 = dot(distance_vector,distance_vector);
			float c = ( INTERACTION_RADIUS - pi.neighbour_distance[i] );
			float pterm = -0.5f * c * W_GRAD_PRESSURE * (pi.pressure + pj.pressure) / pi.neighbour_distance[i];
			float dterm = c * pi.density * pj.density;
			float vterm = W_LAPLACIAN_VISCOSITY * VISCOSITY;
			force += ( pterm * distance_vector + vterm * (pj.velocity_eval - pi.velocity_eval) ) * dterm;
		}
		pi.force = force;
	}
}

void ParticleSystem::moveParticles()
{
	
	for(particle& pi : particles)
	{

		vec2 acceleration = pi.force*PARTICLE_MASS;
		
		//leapfrog integration
		vec2 velocity_next = pi.velocity + acceleration*TIME_STEP;
		pi.velocity_eval = (2.0*pi.velocity + acceleration*TIME_STEP)*0.5;
		pi.velocity = velocity_next;
		velocity_next *= TIME_STEP/(SIM_SCALE);
		pi.position += velocity_next;
	}

}

void ParticleSystem::advance()
{
	updateGrid();
	updateNeighbours();
	calculatePressure();
	calculateSPHForce();
	moveParticles();
	advance_call++;
}
