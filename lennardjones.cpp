#include "lennardjones.h"
#include "system.h"
#include "atom.h"
#include <math.h>
#include "math/vec3.h"

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

double LennardJones::U(Atom i, Atom j)
{
    double rij = distance(i.position,j.position);
    if (rij >= m_sigma) return 0;
    else {
        double elem = this->m_sigma / rij;
        double U = 4 * this->m_epsilon * (pow(elem,12.)-pow(elem,6.));
        return U;
    }
}

void LennardJones::calculateForces(System &system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
}
