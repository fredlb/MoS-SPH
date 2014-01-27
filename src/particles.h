#include <vector>

struct particle
{
	particle(float x, float y, float u, float v)
	{
		m_x = x;
		m_y = y;

		m_u = u;
		m_v = v;
	};

	void update(float dt)
	{
		m_x = m_x + dt*m_u;
		m_y = m_y + dt*m_v;

	};

	float m_x;
	float m_y;

	float m_u; //x-velocity
	float m_v; //y-velocity


};


typedef std::vector<particle> particleSystem;