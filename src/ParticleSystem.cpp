#include "ParticleSystem.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "vec2.h"

#define PI 3.1415926535f

#define MAX_PARTICLES 16384
#define MAX_BORDER_PARTICLES 1000
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


#define STIFFNESS 5.0 //
#define VISCOSITY 0.2 // pascal-seconds
//#define PARTICLE_MASS 0.00020543 //kg


#define REST_DENSITY 600.0 //kg / m^3
#define VEL_LIMIT 200.0 //velocity limit (m/s)
#define PARTICLE_RADIUS 0.004 // m
#define EXT_DAMP 512.0

float PARTICLE_MASS =  0.0f;

const float W_DEFAULT = 315.0f / (64.0f * 3.141592 * pow( INTERACTION_RADIUS, 9) );
const float	W_GRAD_PRESSURE = -45.0f / (3.141592 * pow( INTERACTION_RADIUS, 6) );
const float W_LAPLACIAN_VISCOSITY = 45.0f / (3.141592 * pow( INTERACTION_RADIUS, 6) );

ParticleSystem::ParticleSystem(void)
{
	particles.resize(MAX_PARTICLES);
	grid.resize(GRID_WIDTH*GRID_HEIGHT);
	createParticleField();
	createBorderParticles();
	updateGrid();
	calculateMass();
	std::cout << "Particle mass: " << PARTICLE_MASS << std::endl;
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
	particles.resize(MAX_PARTICLES);

	int rowcolSize = sqrt(MAX_PARTICLES);
	for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
	{
		float stepLength = 0.0038;
		for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
		{
			particles[particleIndexCol + rowcolSize*particleIndexRow].position.x = BORDER_LEFT+0.002 + particleIndexCol*stepLength;
			particles[particleIndexCol + rowcolSize*particleIndexRow].position.y = BORDER_BOTTOM+0.002 + particleIndexRow*stepLength;
			particles[particleIndexCol + rowcolSize*particleIndexRow].is_static = false;
		}
	}
}

void ParticleSystem::createBorderParticles()
{
	border_particles.resize(MAX_BORDER_PARTICLES);
	for( int i = 0; i < MAX_BORDER_PARTICLES; ++i)
	{
		float stepLength = (BORDER_RIGHT - BORDER_LEFT)/MAX_BORDER_PARTICLES;
		border_particles[i].position.x = BORDER_LEFT + i*stepLength;
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
	#pragma omp parallel for
	for(int i = 0; i < MAX_PARTICLES; i++)
	{
		particle& pi = particles[i];
		pi.neighbour_count = 0;
		size_t gi = pi.grid_x;
		size_t gj = pi.grid_y*GRID_WIDTH;

		
		//Loop over border
		/*
		for(int j=0; j<MAX_BORDER_PARTICLES; j++)
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
					std::cout << "I'm on the border" << std::endl;
				}
			}
		}
		*/
		
		//loop over adjacent cells
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
	#pragma omp parallel for
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
				float c = (INTERACTION_RADIUS2 - distance2);
				//#pragma omp atomic
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
	#pragma omp parallel for 
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


	for(particle& pi : particles)
	{
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

	if(particles.size() < MAX_PARTICLES)
	{
		particles.push_back(p);
	}
	else
	{
		particles[draw_counter%particles.size()] = p;
	}
}