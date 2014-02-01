//#pragma once
#include <vector>
#include "vec2.h"

class ParticleSystem
{
public:
	ParticleSystem(void);
	~ParticleSystem(void);
	
	void particlesInit();
	void borderParticlesInit();
	float calculateMass();

	void updateGrid();
	void calculatePressure();
	void calculateForces();
	void advance();

	struct particle
	{
		float m_x;
		float m_y;
		float m_u; //x-velocity
		float m_v; //y-velocity

		float m_massDensity;
		float m_pressure;

		float m_mass;

		particle* next;
	};
	std::vector<particle> particles;
	std::vector<particle> borderParticles;

private:

	#define kMaxNeighbourCount 64
	struct Neighbours
	{
		const particle* particles[kMaxNeighbourCount];
		float r[kMaxNeighbourCount];
		size_t count;
	};
	float kDt;
	bool firstIteration;

	float interactionRadius;
	float IR2;
	float cellSize;
	std::vector<particle*> grid;
	size_t gridWidth;
	size_t gridHeight;
	size_t gridCellCount;

	std::vector<Neighbours> neighbours;
	std::vector<int> gridCoords;

	std::vector<float> vhx;
	std::vector<float> vhy;

	float particleMass;

	float Wdefault(float distance2);
	float kWdefault;
	float kWgradDefault;
	float kWlaplacianDefault;
	float kWlaplacianViscosity;
	float kWgradPressure;

	int averageParticles;

	float viscosityConstant;
	float restDensity;
	float surfaceLimit;
	float surfaceTension;

	float stiffness;

	float g;
};

