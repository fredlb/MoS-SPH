//#pragma once
#include <vector>
#include "vec2.h"

class ParticleSystem
{
public:
	ParticleSystem(void);
	~ParticleSystem(void);

	std::vector<vec2> getParticleCoordinates();
	std::vector<vec2> getParticleCoordinatesNeighbours();
	std::vector<vec2> getParticleCoordinatesGrid();
#define PRESSURE_UNDER 1.0
#define PRESSURE_OVER -1.0
	std::vector<vec2> getParticleCoordinatesPressure(float dir, float limit);
	void advance();
	long int draw_counter;
	void drawParticle(float x, float y, bool is_static);
	void reloadParticleSystem();

private:
	#define MAX_NEIGHBOURS 64
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

		//is moving particle
		bool is_static;
	};
	std::vector<particle> particles;
	std::vector<particle> border_particles;
	std::vector<particle*> grid;

	void createParticleField();
	void updateGrid();
	void updateNeighbours();
	void updateNeighbours2();
	void calculatePressure();
	void calculateSPHForce();
	void calculateMass();
	void moveParticles();
	void createBorderParticles();
	

	long int advance_call;
};

