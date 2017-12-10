#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include "atom.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    for(Atom *atom : this->atoms()) {
        for (int i = 0; i < 3; i++)
        {
            double size = this->m_systemSize[i] / 2;
            if (atom->position[i] < - size / 2) atom->position += size;
            if (atom->position[i] >= size / 2) atom->position -= size;
        }
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum
    double total_momentum[3] ={0.0,0.0,0.0};
    for(Atom *atom : this->atoms()) {
        for (int i = 0; i <3; i++){
            total_momentum[i] = atom->mass()*atom->velocity[i];
        }
    }
    //remove momentum equally on each atom
    for(Atom *atom : this->atoms()){
        double coefficient = this->m_atoms.size()*atom->mass();
        for(int i =0; i<3; i++){
            atom->velocity[i] -= coefficient * total_momentum[i];
        }
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    for( int i = 0; i < numberOfUnitCellsEachDimension; i++){
        for(int j=0; j < numberOfUnitCellsEachDimension; j++){
            for(int k =0; k < numberOfUnitCellsEachDimension; k++){
                vec3 origin(i*latticeConstant,j*latticeConstant,k*latticeConstant);
                Atom *r1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *r2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *r3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *r4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                double b_over_2 = latticeConstant /2.;
                r1->position.set(origin[0],origin[1],origin[2]);
                r2->position.set(origin[0]+b_over_2,origin[1]+b_over_2,origin[2]);
                r3->position.set(origin[0],origin[1]+b_over_2,origin[2]+b_over_2);
                r4->position.set(origin[0]+b_over_2,origin[1],origin[2]+b_over_2);
                r1->resetVelocityMaxwellian(temperature);
                r2->resetVelocityMaxwellian(temperature);
                r3->resetVelocityMaxwellian(temperature);
                r4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(r1);                m_atoms.push_back(r2);
                m_atoms.push_back(r3);                m_atoms.push_back(r4);
            }
        }
    }
    double size_s = numberOfUnitCellsEachDimension * latticeConstant;
    setSystemSize(vec3(size_s, size_s , size_s)); // Remember to set the correct system size!
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
