#include "ParticleSystem.h"
#include <random>
#include <cmath>
#include <iostream>
#include <vector>
#include "vec2.h"

#define kPi 3.14159265359


ParticleSystem::ParticleSystem(void)
{

	kDt = 0.00025f;
	firstIteration = true;
	
	interactionRadius = 0.035f;
	IR2 = interactionRadius*interactionRadius;
	cellSize = (2.0f*interactionRadius);
	
	gridWidth = (size_t)(2.0f / cellSize);
	gridHeight = (size_t)(2.0f / cellSize);
	gridCellCount = gridWidth * gridHeight;
	
	grid.resize(gridCellCount);
	
	particlesInit();
	borderParticlesInit();

	gridCoords.resize(particles.size()*2);
	neighbours.resize(particles.size());
	vhx.resize(particles.size());
	vhy.resize(particles.size());

	
	updateGrid();
	
	kWdefault = (315/(64*kPi*pow(interactionRadius,9)));
	kWgradDefault = -(945/(32*kPi*pow(interactionRadius,9)));
	kWlaplacianDefault = -(945/(32*kPi*pow(interactionRadius,9)));
	kWlaplacianViscosity = (45/(kPi*pow(interactionRadius,6)));
	kWgradPressure = -(45/(kPi*pow(interactionRadius,6)));

	averageParticles = 20;

	viscosityConstant = 3.5f;
	restDensity = 988.0f;
	surfaceLimit = sqrt(restDensity/averageParticles);
	surfaceTension = 0.0728f;

	stiffness = 100.0f;

	particleMass = calculateMass();

	g = -9.81;
	
}


ParticleSystem::~ParticleSystem(void)
{
}

std::vector<vec2> ParticleSystem::getCoordinateVector()
{
	std::vector<vec2> coordinateVector(particles.size()+borderParticles.size());
	for(int i = 0; i < particles.size(); ++i)
	{
		vec2 p;
		p.x = particles[i].m_x;
		p.y = particles[i].m_y;
		coordinateVector[i] = p;
	}
	//#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < borderParticles.size(); ++i)
	{
		vec2 p;
		p.x = borderParticles[i].m_x;
		p.y = borderParticles[i].m_y;
		coordinateVector[particles.size() + i] = p;
	}
	return coordinateVector;
}

#define MAX_PARTICLES 1024
#define kOffset 0.5f
void ParticleSystem::particlesInit()
{
	std::mt19937 eng((std::random_device())());
	std::uniform_real_distribution<> pos_dist(-0.001,0.001);
	particles.resize(MAX_PARTICLES);
	int rowcolSize = sqrt(particles.size());

	for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
	{
		float stepLength = 1.0f/rowcolSize;
		for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
		{
			particles[particleIndexCol + rowcolSize*particleIndexRow].m_x = -kOffset + particleIndexCol*stepLength; // + pos_dist(eng);
			particles[particleIndexCol + rowcolSize*particleIndexRow].m_y = -kOffset + particleIndexRow*stepLength; // + pos_dist(eng);
		}
	}
}

#define MAX_BORDER 300
void ParticleSystem::borderParticlesInit()
{
	borderParticles.resize(MAX_BORDER);
	float stepLengthx = 2.0f/borderParticles.size();

	for(int i = 0; i < borderParticles.size(); i++)
	{
		borderParticles[i].m_x = -1.0f + stepLengthx*i;
		borderParticles[i].m_y = -0.98f;
		borderParticles[i].m_mass = 1;
		borderParticles[i].m_massDensity = 10*restDensity;
		borderParticles[i].m_pressure = stiffness*(9*restDensity);
	}
}

float ParticleSystem::calculateMass()
{
	float density = 0.0f; 


	for(int i = 0; i < particles.size(); ++i)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;

		size_t gi = gridCoords[i*2];
		size_t gj = gridCoords[i*2+1]*gridWidth;

		//loop over cells 
		for (size_t ni=gi-1; ni<=gi+1; ++ni)
		{
			for (size_t nj=gj-gridWidth; nj<=gj+gridWidth; nj+=gridWidth)
			{
				//loop over neighbors
				for (particle* ppj=grid[ni+nj]; NULL!=ppj; ppj=ppj->next)
				{
					const particle& pj = *ppj;
					//std::cout << "pi x:" << pi.m_x << std::endl;

					//std::cout << "pj x:" << pj.m_x << std::endl;
					float dx = pi.m_x - pj.m_x;
					//std::cout << "dx: " << dx << std::endl;
					float dy = pi.m_y - pj.m_y;
					float distance2 = dx*dx + dy*dy;
					if(distance2 < interactionRadius*interactionRadius)
					{
						density += kWdefault * (interactionRadius*interactionRadius - distance2)*(interactionRadius*interactionRadius - distance2)*(interactionRadius*interactionRadius - distance2);
					}
				}
			}
		}
	}
	float dA = density/particles.size();
	float mass = (dA*restDensity)/(dA*dA);
	for(size_t i = 0; i < particles.size(); ++i)
	{
		particle& pi = particles[i];
		pi.m_mass = mass;
	}
	return mass;
}



void ParticleSystem::updateGrid()
{
	memset(&grid[0], 0, grid.size()*sizeof(grid[0]));
	
	for(size_t i = 0; i < particles.size(); i++)
	{
		particle& pi = particles[i];
		
		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;
		
		if (x < 1)
			x = 1;
		else if (x > gridWidth-2)
			x = gridWidth-2;
		
		if (y < 1)
			y = 1;
		else if (y > gridHeight-2)
			y = gridHeight-2;
		
		pi.next = grid[x+y*gridWidth];
		grid[x+y*gridWidth] = &pi;
		
		gridCoords[i*2] = x;
		gridCoords[i*2+1] = y;
	}
}

void ParticleSystem::calculatePressure()
{
	//Mass-density and pressure loop


	for(int i = 0; i < particles.size(); ++i)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;

		size_t gi = gridCoords[i*2];
		size_t gj = gridCoords[i*2+1]*gridWidth;

		float massDensity = 0.0f;
		neighbours[i].count = 0;

		//Loop over border
		for(size_t j = 0; j < borderParticles.size(); j++)
		{
			particle& bp = borderParticles[j];
			float pm = bp.m_mass;

			float dx = pi.m_x - bp.m_x;
			float dy = pi.m_y - bp.m_y;
			float distance2 = dx*dx + dy*dy;
			if(distance2 < IR2)
			{
				//Density
				massDensity += particleMass*Wdefault(distance2);
				if(neighbours[i].count < kMaxNeighbourCount)
				{
					neighbours[i].particles[neighbours[i].count] = &bp;
					neighbours[i].r[neighbours[i].count] = sqrt(distance2);
					++neighbours[i].count;
					//std::cout << "I'm on the border" << std::endl;
				}
			}
		}
		//loop over cells
		for (int ni=gi-1; ni<=gi+1; ++ni)
		{
			for (int nj=gj-gridWidth; nj<=gj+gridWidth; nj+=gridWidth)
			{
				//loop over neighbors
				for (particle* ppj=grid[ni+nj]; NULL!=ppj; ppj=ppj->next)
				{
					//do fancy math
					//std::cout << "ppj x: " << ppj->m_x << std::endl;
					float dx = pi.m_x - ppj->m_x;
					float dy = pi.m_y - ppj->m_y;
					float distance2 = dx*dx + dy*dy;

					if(distance2 < IR2)
					{
						//Density
						massDensity += particleMass*kWdefault* (IR2 - distance2)*(IR2 - distance2)*(IR2 - distance2);
						//std::cout << particleMass;
						if(neighbours[i].count < kMaxNeighbourCount)
						{
							neighbours[i].particles[neighbours[i].count] = ppj;
							neighbours[i].r[neighbours[i].count] = sqrt(distance2);
							++neighbours[i].count;
						}
					}
				}
			}
		}
		//save massDensity
		pi.m_massDensity = massDensity;
		//std::cout << massDensity << std::endl;
		//save Pressure
		pi.m_pressure = stiffness * (massDensity - restDensity);
		//std::cout << pi.m_pressure << std::endl;
	}



}

void ParticleSystem::calculateForces()
{
	//Force loop
	//
	for(int i = 0; i < particles.size(); ++i)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;

		float pressureForcex = 0.0f; 
		float pressureForcey = 0.0f;
		float viscosityForcex = 0.0f;
		float viscosityForcey = 0.0f;
		float normalx = 0.0f;
		float normaly = 0.0f;
		float gradNormal = 0.0f;
		float surfaceTensionForcex = 0.0f;
		float surfaceTensionForcey = 0.0f;

		size_t gi = gridCoords[i*2];
		size_t gj = gridCoords[i*2+1]*gridWidth;

		float mdi = pi.m_massDensity;
		float mdj = 0.0f;
		float gravity = g*mdi;
		//loop over neighbours 
				for(int j=0; j < neighbours[i].count; ++j)
				{
					const particle* ppj = neighbours[i].particles[j];
					float massi = pi.m_mass;
					float massj = ppj->m_mass;
					float dx = pi.m_x - ppj->m_x;
					float dy = pi.m_y - ppj->m_y;
					float distance2 = dx*dx + dy*dy;
					if(distance2 < IR2)
					{
						mdj = ppj->m_massDensity;
						
						normalx += (particleMass/mdj)*kWgradDefault*(IR2-distance2)*(IR2-distance2)*dx;
						normaly += (particleMass/mdj)*kWgradDefault*(IR2-distance2)*(IR2-distance2)*dy;

						gradNormal += (particleMass/mdj)*kWlaplacianDefault*(IR2-distance2)*(IR2-distance2)*(3*IR2-7*distance2);

						if( distance2 != 0)
						{
							float distance = sqrt(distance2);
							// uj - ui
							float velocityDiffu = ppj->m_u - pi.m_u;
							float velocityDiffv = ppj->m_v - pi.m_v;

							float pressi = pi.m_pressure;
							float pressj = ppj->m_pressure;
							pressureForcex += ((pressi/pow(mdi,2))+(pressj/pow(mdj,2)))*particleMass*(interactionRadius-distance)*(interactionRadius-distance)*(dx/distance)*kWgradPressure;
							pressureForcey += ((pressi/pow(mdi,2))+(pressj/pow(mdj,2)))*particleMass*(interactionRadius-distance)*(interactionRadius-distance)*(dy/distance)*kWgradPressure;
							
							viscosityForcex += velocityDiffu * (particleMass/mdj) * kWlaplacianViscosity*(interactionRadius-distance);
							viscosityForcey += velocityDiffv * (particleMass/mdj) * kWlaplacianViscosity*(interactionRadius-distance);

						}
					}
				}
		float normalLenght = sqrt(normalx*normalx + normaly*normaly);
		if(normalLenght > surfaceLimit){
			surfaceTensionForcex = - surfaceTension  * gradNormal * (normalx/normalLenght);
			surfaceTensionForcey = - surfaceTension  * gradNormal * (normaly/normalLenght);
		}
		pressureForcex = -mdi * pressureForcex;
		pressureForcey = -mdi * pressureForcey;
		
		viscosityForcex = viscosityConstant*viscosityForcex;
		viscosityForcey = viscosityConstant*viscosityForcey;

		float accelerationX = (pressureForcex + viscosityForcex + surfaceTensionForcex)/mdi;
		float accelerationY = (pressureForcey + viscosityForcey + surfaceTensionForcey + gravity)/mdi;

		//True leap-frog
                if(firstIteration){
                        vhx[i] = pi.m_u + 0.5*accelerationX*kDt;
                        vhy[i] = pi.m_v + 0.5*accelerationY*kDt;
                        
                        pi.m_u += accelerationX*kDt;
                        pi.m_v += accelerationY*kDt;

                        pi.m_x += vhx[i]*kDt;
                        pi.m_y += vhy[i]*kDt;
						firstIteration = false;
                }else{
                        vhx[i] += accelerationX*kDt;
                        vhy[i] += accelerationY*kDt;
                        
                        pi.m_u = vhx[i] + 0.5*accelerationX*kDt;
                        pi.m_v = vhy[i] + 0.5*accelerationY*kDt;
                        pi.m_x += vhx[i]*kDt;
                        pi.m_y += vhy[i]*kDt;

                }

	}


}

void ParticleSystem::advance()
{
	updateGrid();
	calculatePressure();
	calculateForces();
}

float ParticleSystem::Wdefault(float distance2)
{
	float W = (315/(64*kPi*pow(interactionRadius,9))) * pow((pow(interactionRadius,2) -distance2),3);
	return W;
}