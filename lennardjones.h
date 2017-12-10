#ifndef LENNARDJONES_H
#define LENNARDJONES_H
#include "atom.h"

class LennardJones
{
private:
    double m_sigma = 1.0;
    double m_epsilon = 1.0;
    double m_potentialEnergy = 0;

public:
    LennardJones() { }
    void calculateForces(class System &system);
    double U(Atom i, Atom j);
    double potentialEnergy() const;
    double sigma() const;
    void setSigma(double sigma);
    double epsilon() const;
    void setEpsilon(double epsilon);
};
#endif
