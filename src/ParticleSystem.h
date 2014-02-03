//#pragma once
#include <vector>
#include "vec2.h"

using namespace std;

class ParticleSystem
{
public:
	ParticleSystem(void);
	~ParticleSystem(void);

	std::vector<vec2> getParticleCoordinates();
	void advance();
	void moveParticleTo(float x, float y);
	void moveParticleTo(vec2 xy);

private:
	#define MAX_NEIGHBOURS 128
	struct particle
	{
		//physical properties
		vec2 position;
		vec2 velocity;
		vec2 velocity_eval;
		vec2 force;
		float pressure;
		float density;

		//neighbour list
		particle* neighbours[MAX_NEIGHBOURS];
		float neighbour_distance[MAX_NEIGHBOURS];
		size_t neighbour_count;


		//for grid linked list
		particle* next;
		//grid coordinates
		size_t grid_x;
		size_t grid_y;
	};
	vector<particle> particles;
	vector<particle> border_particles;
	vector<particle*> grid;

	void createParticleField();
	void updateGrid();
	void updateNeighbours();
	void calculatePressure();
	void calculateSPHForce();
	void moveParticles();

	long int advance_call;
};

