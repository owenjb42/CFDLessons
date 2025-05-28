#include "BoundaryControl.hpp"
#include "Solver.hpp"

/////////////////
// Fixed Inlet //
/////////////////

void FixedInletBoundaryCondition::ApplyForVelocity(Solver& solver)
{
    int offset = direction == 0 ? 1 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        double flux = inlet_velocity * solver.dy;
        double coeff = solver.fluid.density * flux / (solver.dx * solver.dy);
        for (auto [i, j] : boundary_faces)
        {
            solver.p_scr(i + cell_offset, j) += -inlet_velocity * solver.dy;
            solver.u_coeff(i + offset, j) += coeff;
            solver.u_scr(i + offset, j) += coeff * offset * inlet_velocity;
        }
    }
    else // v
    {
        double flux = inlet_velocity * solver.dx;
        double coeff = solver.fluid.density * flux / (solver.dx * solver.dy);
        for (auto [i, j] : boundary_faces)
        {
            solver.p_scr(i, j + cell_offset) += -inlet_velocity * solver.dx;
            solver.v_coeff(i, j + offset) += coeff;
            solver.v_scr(i, j + offset) += coeff * offset * inlet_velocity;
        }
    }
}

void FixedInletBoundaryCondition::ApplyForTemperature(Solver& solver)
{
    int offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        double massflux = inlet_velocity * solver.fluid.density * solver.dy;
        double coeff = massflux / (solver.dx * solver.dy);

        for (auto [i, j] : boundary_faces)
        {
            solver.t_coeff(i + offset, j) += coeff;
            solver.t_scr(i + offset, j) += coeff * inlet_temperature;
        }
    }
    else // v
    {
        double massflux = inlet_velocity * solver.fluid.density * solver.dx;
        double coeff = massflux / (solver.dx * solver.dy);

        for (auto [i, j] : boundary_faces)
        {
            solver.t_coeff(i, j + offset) += coeff;
            solver.t_scr(i, j + offset) += coeff * inlet_temperature;
        }
    }
}

void FixedInletBoundaryCondition::CorrectBoundaryCellVelocities(Solver& solver)
{
    int offset = direction == 0 ? 1 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            solver.cell_u(i + cell_offset, j) += offset * inlet_velocity / 2.0;
        }
    }
    else // v
    {
        int cell_offset = direction == 0 ? 0 : -1;
        for (auto [i, j] : boundary_faces)
        {
            solver.cell_v(i, j + cell_offset) += offset * inlet_velocity / 2.0;
        }
    }
}

//////////////////
// Fixed Outlet //
//////////////////

void FixedOutletBoundaryCondition::ApplyForVelocity(Solver& solver)
{
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            int offset = direction == 0 ? -1 : 1;
            int cell_offset = direction == 0 ? 0 : -1;
            solver.p_scr(i + cell_offset, j) += outlet_velocity * solver.dy;
        }
    }
    else // v
    {
        for (auto [i, j] : boundary_faces)
        {
            int offset = direction == 0 ? -1 : 1;
            int cell_offset = direction == 0 ? 0 : -1;
            solver.p_scr(i, j + cell_offset) += outlet_velocity * solver.dx;
        }
    }
}
void FixedOutletBoundaryCondition::CorrectBoundaryCellVelocities(Solver& solver)
{
    int offset = direction == 0 ? -1 : 1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            solver.cell_u(i + cell_offset, j) += offset * outlet_velocity / 2.0;
        }
    }
    else // v
    {
        int cell_offset = direction == 0 ? 0 : -1;
        for (auto [i, j] : boundary_faces)
        {
            solver.cell_v(i, j + cell_offset) += offset * outlet_velocity / 2.0;
        }
    }
}

///////////////////
// Open Boundary //
///////////////////
void OpenBoundaryCondition::ApplyForVelocity(Solver& solver)
{
    double fix_coef = 1000.0;

    int offset = direction == 0 ? 0 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            double coeff = fix_coef * solver.dy;
            double dp = (pressure - solver.p(i + cell_offset, j));
            double u_b = dp / solver.fluid.density;
            solver.p_scr(i + cell_offset, j) += coeff * -u_b;
            solver.p_coef(i + cell_offset, j) += coeff;

            double massflux = fix_coef * (pressure - solver.p(i + cell_offset, j)) * solver.dy;
            if (massflux > 0.0)
            {
                double coeff = massflux / (solver.dx * solver.dy);
                solver.u_coeff(i + offset, j) += coeff;
                solver.u_scr(i + offset, j) += coeff * offset * u_b;
            }
        }
    }
    else // v
    {
        for (auto [i, j] : boundary_faces)
        {
            double coeff = fix_coef * solver.dx;
            double dp = (pressure - solver.p(i, j + cell_offset));
            double u_b = dp / solver.fluid.density;
            solver.p_scr(i, j + cell_offset) += coeff * -u_b;
            solver.p_coef(i, j + cell_offset) += coeff;

            double massflux = fix_coef * (pressure - solver.p(i, j + cell_offset)) * solver.dx;
            if (massflux > 0.0)
            {
                double coeff = massflux / (solver.dx * solver.dy);
                solver.v_coeff(i, j + offset) += coeff;
                solver.v_scr(i, j + offset) += coeff * offset * u_b;
            }
        }
    }
}
void OpenBoundaryCondition::ApplyForTemperature(Solver& solver)
{
    double fix_coef = 1000.0;

    int offset = direction == 0 ? 0 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            double massflux = fix_coef * (pressure - solver.p(i + cell_offset, j)) * solver.dy;
            if (massflux > 0.0)
            {
                double coeff = massflux / (solver.dx * solver.dy);
                solver.t_coeff(i + offset, j) += coeff;
                solver.t_scr(i + offset, j) += coeff * temperature;
            }
        }
    }
    else // v
    {
        for (auto [i, j] : boundary_faces)
        {
            double massflux = fix_coef * (pressure - solver.p(i, j + cell_offset)) * solver.dx;
            if (massflux > 0.0)
            {
                double coeff = massflux / (solver.dx * solver.dy);
                solver.t_coeff(i, j + offset) += coeff;
                solver.t_scr(i, j + offset) += coeff * temperature;
            }
        }
    }
}
void OpenBoundaryCondition::CorrectBoundaryCellVelocities(Solver& solver)
{
    double fix_coef = 1000.0;

    int offset = direction == 0 ? 1 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            double dp = (pressure - solver.p(i + cell_offset, j));
            double u_b = fix_coef * dp / solver.fluid.density;
            solver.cell_u(i + cell_offset, j) += offset * u_b / 2.0;
        }
    }
    else // v
    {
        for (auto [i, j] : boundary_faces)
        {
            double dp = (pressure - solver.p(i, j + cell_offset));
            double u_b = fix_coef * dp / solver.fluid.density;
            solver.cell_v(i, j + cell_offset) += offset * u_b / 2.0;
        }
    }
}

///////////////////////
// Friction Boundary //
///////////////////////
void FrictionBoundaryCondition::ApplyFriction(Solver& solver)
{

}
