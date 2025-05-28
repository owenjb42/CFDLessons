#pragma once
#include <vector>

class Solver;

struct BoundaryCondition
{
    BoundaryCondition(bool direction, bool component) : direction(direction), component(component) {}

    bool direction; // 0 : +, 1 : -
    bool component; // 0 : u, 1 : v
    std::vector<std::pair<int, int>> boundary_faces;
};

struct FixedInletBoundaryCondition : public BoundaryCondition
{
    FixedInletBoundaryCondition(bool direction, bool component) : BoundaryCondition(direction, component) {}

    void ApplyForVelocity(Solver& solver);
    void ApplyForTemperature(Solver& solver);
    void CorrectBoundaryCellVelocities(Solver& solver);

    double inlet_temperature{};
    double inlet_velocity{};
};

struct FixedOutletBoundaryCondition : public BoundaryCondition
{
    FixedOutletBoundaryCondition(bool direction, bool component) : BoundaryCondition(direction, component) {}

    void ApplyForVelocity(Solver& solver);
    void CorrectBoundaryCellVelocities(Solver& solver);

    double outlet_velocity{};
};

struct OpenBoundaryCondition : public BoundaryCondition
{
    OpenBoundaryCondition(bool direction, bool component) : BoundaryCondition(direction, component) {}

    void ApplyForVelocity(Solver& solver);
    void ApplyForTemperature(Solver& solver);
    void CorrectBoundaryCellVelocities(Solver& solver);

    double temperature{};
    double pressure{};
};

struct FrictionBoundaryCondition
{
    void ApplyFriction(Solver& solver);
    double wall_velocity = 0.0;
};