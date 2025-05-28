#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cmath>
#include <limits>
#include <thread>
#include <algorithm>

#include "Fields/PhysicalField.hpp"
#include "Fields/FlagField.hpp"
#include "GUI/Interface.hpp"
#include "BoundaryControl.hpp"

class Solver
{
public:
    Solver(Interface& interface) : interface(interface), nx(interface.nx), ny(interface.ny), dx(interface.dx), dy(interface.dy)
    {
        SetBlockedFaces(interface);
        SetBoundaries(interface);
        interface.SetData(*this);
    }
    bool is_solving{ true };
    Interface& interface;

    void ComputeDivergence()
    {
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                divergence(i, j) = ((u(i + 1, j) - u(i, j)) * dy +
                    (v(i, j + 1) - v(i, j)) * dx +
                    p_scr(i, j));
            }
        }
    }

    void ComputeConvectiveTerms()
    {
        double flux = 0.0;
        double value = 0.0;
        double coeff = 0.0;

        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (!u_face_flags(i, j, Flag::Open))
                    continue;

                if (u_face_flags(i - 1, j, Flag::Open))
                {
                    // Inline Left: u flux
                    flux = dy * ((u(i, j) + u(i - 1, j)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i - 1, j) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (u_face_flags(i + 1, j, Flag::Open))
                {
                    // Inline Right: u flux
                    flux = dy * -((u(i, j) + u(i + 1, j)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i + 1, j) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (j != 0 && v_face_flags(i, j, Flag::Open) && v_face_flags(i - 1, j, Flag::Open) && u_face_flags(i, j - 1, Flag::Open))
                {
                    // Bottom Left: u flux
                    flux = (dx / 2.0) * v(i - 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: u flux
                    flux = (dx / 2.0) * v(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (j != ny - 1 && v_face_flags(i, j + 1, Flag::Open) && v_face_flags(i - 1, j + 1, Flag::Open) && u_face_flags(i, j + 1, Flag::Open))
                {
                    // Top Left: u flux
                    flux = -(dx / 2.0) * v(i - 1, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Top Right: u flux
                    flux = -(dx / 2.0) * v(i, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if (!v_face_flags(i, j, Flag::Open))
                    continue;

                if (v_face_flags(i, j - 1, Flag::Open))
                {
                    // Inline Left: v flux
                    flux = dy * ((v(i, j) + v(i, j - 1)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j - 1) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (v_face_flags(i, j + 1, Flag::Open))
                {
                    // Inline Right: v flux
                    flux = dy * -((v(i, j) + v(i, j + 1)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j + 1) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (i != 0 && u_face_flags(i + 1, j, Flag::Open) && u_face_flags(i + 1, j - 1, Flag::Open) && v_face_flags(i + 1, j, Flag::Open))
                {
                    // Bottom Left: v flux
                    flux = (dy / 2.0) * u(i, j - 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: v flux
                    flux = (dy / 2.0) * u(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (i != nx - 1 && u_face_flags(i, j, Flag::Open) && u_face_flags(i, j - 1, Flag::Open) && v_face_flags(i - 1, j, Flag::Open))
                {
                    // Top Left: v flux
                    flux = -(dy / 2.0) * u(i + 1, j - 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Top Right: v flux
                    flux = -(dy / 2.0) * u(i + 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }
            }
        }
    }

    void ComputeDiffusiveTerms()
    {
        double coeff = 0.0;
        double value = 0.0;
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (!u_face_flags(i, j, Flag::Open))
                    continue;

                if (u_face_flags(i - 1, j, Flag::Open))
                {
                    // Inline Left: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i - 1, j);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (u_face_flags(i + 1, j, Flag::Open))
                {
                    // Inline Right: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i + 1, j);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != 0 && v_face_flags(i, j, Flag::Open) && v_face_flags(i - 1, j, Flag::Open) && u_face_flags(i, j - 1, Flag::Open))
                {
                    // Bottom: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = u(i, j - 1);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != ny - 1 && v_face_flags(i, j + 1, Flag::Open) && v_face_flags(i - 1, j + 1, Flag::Open) && u_face_flags(i, j + 1, Flag::Open))
                {
                    // Top: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = u(i, j + 1);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if (!v_face_flags(i, j, Flag::Open))
                    continue;

                if (v_face_flags(i, j - 1, Flag::Open))
                {
                    // Inline Left: v diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = v(i, j - 1);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (v_face_flags(i, j + 1, Flag::Open))
                {
                    // Inline Right: v diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = v(i, j + 1);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (i != 0 && u_face_flags(i, j, Flag::Open) && u_face_flags(i, j - 1, Flag::Open) && v_face_flags(i - 1, j, Flag::Open))
                {
                    // Bottom: v diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = v(i - 1, j);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (i != nx - 1 && u_face_flags(i + 1, j, Flag::Open) && u_face_flags(i + 1, j - 1, Flag::Open) && v_face_flags(i + 1, j, Flag::Open))
                {
                    // Top: v diff
                    coeff = fluid.kinematic_viscosity/ (dx * dx);
                    value = v(i + 1, j);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }
            }
        }
    }

    void ComputeConvectiveFluxesAndDiffusiveTermsForTemperature()
    {
        double value = 0.0;
        double coeff = 0.0;
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (u_face_flags(i, j, Flag::Open))
                {
                    value = t(i - 1, j);

                    // Left: Diff
                    coeff = (fluid.conductivity / fluid.cp) / (dx * dx);
                    t_coeff(i, j) += coeff;
                    t_scr(i, j) += coeff * value;

                    // Left: Cnv
                    coeff = fluid.density * dy * u(i, j) / (dy * dx);
                    if (coeff > 0.0)
                    {
                        t_coeff(i, j) += coeff;
                        t_scr(i, j) += coeff * value;
                    }
                }

                if (u_face_flags(i + 1, j, Flag::Open))
                {
                    value = t(i + 1, j);

                    // Right: Diff
                    coeff = (fluid.conductivity / fluid.cp) / (dx * dx);
                    t_coeff(i, j) += coeff;
                    t_scr(i, j) += coeff * value;

                    // Right: Cnv
                    coeff = -fluid.density * dy * u(i + 1, j) / (dy * dx);
                    if (coeff > 0.0)
                    {
                        t_coeff(i, j) += coeff;
                        t_scr(i, j) += coeff * value;
                    }
                }

                if (v_face_flags(i, j, Flag::Open))
                {
                    value = t(i, j - 1);

                    // Bottom: Diff
                    coeff = (fluid.conductivity / fluid.cp) / (dy * dy);
                    t_coeff(i, j) += coeff;
                    t_scr(i, j) += coeff * value;

                    // Bottom: Cnv
                    coeff = fluid.density * dy * v(i, j) / (dy * dx);
                    if (coeff > 0.0)
                    {
                        t_coeff(i, j) += coeff;
                        t_scr(i, j) += coeff * value;
                    }
                }

                if (v_face_flags(i, j + 1, Flag::Open))
                {
                    value = t(i, j + 1);

                    // Top: Diff
                    coeff = (fluid.conductivity / fluid.cp) / (dy * dy);
                    t_coeff(i, j) += coeff;
                    t_scr(i, j) += coeff * value;

                    // Top: Cnv
                    coeff = -fluid.density * dy * v(i, j + 1) / (dy * dx);
                    if (coeff > 0.0)
                    {
                        t_coeff(i, j) += coeff;
                        t_scr(i, j) += coeff * value;
                    }
                }
            }
        }
    }

    void CalculateCellVelocities()
    {
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                cell_u(i, j) = (u(i, j) + u(i + 1, j)) / 2.0;
                cell_v(i, j) = (v(i, j) + v(i, j + 1)) / 2.0;
            }
        }

        CorrectBoundaryCellVelocities();
    }

    void SolvePressure()
    {
        for (int itter = 0; itter < innerPressureItterations; ++itter)
        {
            std::swap(p_correction_old, p_correction);
            double max_correction = 0.0;
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < ny; ++j)
                {
                    double pc = 0.0;

                    double Ae = (u_coeff(i, j) != 0.0) ? (dy / dx) / u_coeff(i, j) : 0.0;
                    double Aw = (u_coeff(i + 1, j) != 0.0) ? (dy / dx) / u_coeff(i + 1, j) : 0.0;
                    double As = (v_coeff(i, j) != 0.0) ? (dx / dy) / v_coeff(i, j) : 0.0;
                    double An = (v_coeff(i, j + 1) != 0.0) ? (dx / dy) / v_coeff(i, j + 1) : 0.0;

                    double Ap = Ae + Aw + As + An + p_coef(i, j);

                    if (Ap == 0) { continue; }

                    if (u_face_flags(i, j, Flag::Open))
                        pc += (Ae / Ap) * p_correction_old(i - 1, j);
                    if (u_face_flags(i + 1, j, Flag::Open))
                        pc += (Aw / Ap) * p_correction_old(i + 1, j);
                    if (v_face_flags(i, j, Flag::Open))
                        pc += (As / Ap) * p_correction_old(i, j - 1);
                    if (v_face_flags(i, j + 1, Flag::Open))
                        pc += (An / Ap) * p_correction_old(i, j + 1);

                    pc -= divergence(i, j) / Ap;

                    p_correction(i, j) = pc;
                    max_correction = std::max(max_correction, std::abs(pc));
                }
            }
            if (relaxation * max_correction < maxPressureResidual) { break; }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                p(i, j) += relaxation * p_correction(i, j);
            }
        }
    }

    void ApplyPressureCorrectionToVelocity()
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (u_coeff(i, j) != 0) { u(i, j) += (p_correction(i - 1, j) - p_correction(i, j)) / (dx * u_coeff(i, j)); }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if (v_coeff(i, j) != 0) { v(i, j) += (p_correction(i, j - 1) - p_correction(i, j)) / (dy * v_coeff(i, j)); }
            }
        }
    }

    void SolveVelocitiesForMomentumEquation()
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (u_coeff(i, j) != 0.0) u(i, j) += relaxation * (u_scr(i, j) - u(i, j) * u_coeff(i, j) + (p(i - 1, j) - p(i, j)) / dx) / (u_coeff(i, j) + fluid.density / dt);
            }
        }
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if (v_coeff(i, j) != 0.0) v(i, j) += relaxation * (v_scr(i, j) - v(i, j) * v_coeff(i, j) + (p(i, j - 1) - p(i, j)) / dy) / (v_coeff(i, j) + fluid.density / dt);
            }
        }
    }

    void SolveTemperatures()
    {
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                t(i, j) = t_scr(i, j) / (t_coeff(i, j) + 10e-20);
            }
        }
    }

    ///////////////////////////////
    // Main time stepping method //
    ///////////////////////////////
    void solve(int num_iterations) 
    {
        is_solving = true;

        for (int iter = 0; iter < num_iterations; ++iter) 
        {
            if (!is_solving)
                break;

            v_scr.reset();
            u_scr.reset();
            p_correction.reset();
            p_scr.reset();
            p_coef.reset();
            u_coeff.reset();
            v_coeff.reset();
            t_scr.reset();
            t_coeff.reset();
            t_previous = t;
            
            // Solve Temperature

            ApplyTemperatureBoundaryConditions();
            
            ComputeConvectiveFluxesAndDiffusiveTermsForTemperature();
            
            SolveTemperatures();

            // Solve Pressure/Velocity

            ApplyVelocityBoundaryConditions();

            ComputeConvectiveTerms();

            ComputeDiffusiveTerms();

            SolveVelocitiesForMomentumEquation();

            ApplyPressureBoundaryConditions();

            ComputeDivergence();

            SolvePressure();

            ApplyPressureCorrectionToVelocity();

            // Check Residuals

            /* ADD residual code here */

            CalculateCellVelocities();
            interface.SetData(*this);
        }

        is_solving = false;
    }

    ////////////////
    // Boundaries //
    ////////////////

    void ApplyVelocityBoundaryConditions() 
    {
        for (auto& boundary : inlet_boundary_conditions)
            boundary.ApplyForVelocity(*this);
        
        for (auto& boundary : outlet_boundary_conditions)
            boundary.ApplyForVelocity(*this);

        for (auto& boundary : open_boundary_conditions)
            boundary.ApplyForVelocity(*this);

        friction_boundary_condition.ApplyFriction(*this);
    }

    void ApplyTemperatureBoundaryConditions()
    {
        for (auto& boundary : inlet_boundary_conditions)
            boundary.ApplyForTemperature(*this);

        for (auto& boundary : open_boundary_conditions)
            boundary.ApplyForTemperature(*this);
    }

    void CorrectBoundaryCellVelocities()
    {
        for (auto& boundary : inlet_boundary_conditions)
            boundary.CorrectBoundaryCellVelocities(*this);

        for (auto& boundary : outlet_boundary_conditions)
            boundary.CorrectBoundaryCellVelocities(*this);

        for (auto& boundary : open_boundary_conditions)
            boundary.CorrectBoundaryCellVelocities(*this);
    }

    void SetBlockedFaces(Interface& interface)
    {
        u_face_flags.setFlag(Flag::Open);
        v_face_flags.setFlag(Flag::Open);

        // Boundary faces
        for (int j = 0; j < ny; ++j)
        {
            u_face_flags.clearFlag(0, j, Flag::Open);
            u_face_flags.clearFlag(nx, j, Flag::Open);

            u_face_flags.setFlag(0, j, Flag::FrictionFace);
            u_face_flags.setFlag(nx, j, Flag::FrictionFace);
        }
        for (int i = 0; i < nx; ++i)
        {
            v_face_flags.clearFlag(i, 0, Flag::Open);
            v_face_flags.clearFlag(i, ny, Flag::Open);

            v_face_flags.setFlag(i, 0, Flag::FrictionFace);
            v_face_flags.setFlag(i, ny, Flag::FrictionFace);
        }

        // User set faces
        for (auto& face : interface.blocked_faces)
        {
            if (face.component_dir == 0)
            {
                u_face_flags.clearFlag(face.i, face.j, Flag::Open);
                u_face_flags.setFlag(face.i, face.j, Flag::FrictionFace);
            }
            else
            {
                v_face_flags.clearFlag(face.i, face.j, Flag::Open);
                v_face_flags.setFlag(face.i, face.j, Flag::FrictionFace);
            }
        }
    }

    void SetBoundaries(Interface& interface)
    {
        // User set 
        for (auto& boundary : interface.boundary_faces)
        {
            if (boundary.component_dir == 0)
            {
                u_face_flags.clearFlag(boundary.i, boundary.j, Flag::FrictionFace);
                u_face_flags.clearFlag(boundary.i, boundary.j, Flag::FrictionFace);
            }
            else
            {
                v_face_flags.clearFlag(boundary.i, boundary.j, Flag::FrictionFace);
                v_face_flags.clearFlag(boundary.i, boundary.j, Flag::FrictionFace);
            }

            if (boundary.type == Boundary::FixedInflow)
            {
                auto& bc = inlet_boundary_conditions.emplace_back(boundary.boundary_dir, boundary.component_dir);
                bc.inlet_temperature = boundary.optional_temp;
                bc.inlet_velocity = boundary.optional_velocity;
                bc.boundary_faces.push_back({ boundary.i, boundary.j });
                boundary.component_dir == 0 ? u_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open) : v_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open);
            }
            if (boundary.type == Boundary::FixedOutflow)
            {
                auto& bc = outlet_boundary_conditions.emplace_back(boundary.boundary_dir, boundary.component_dir);
                bc.outlet_velocity = boundary.optional_velocity;
                bc.boundary_faces.push_back({ boundary.i, boundary.j });
                boundary.component_dir == 0 ? u_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open) : v_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open);
            }
            else if (boundary.type == Boundary::Open)
            {
                if ((boundary.component_dir == 0 && boundary.i != nx) || (boundary.component_dir == 1 && boundary.j != ny))
                {
                    auto& bc = open_boundary_conditions.emplace_back(0, boundary.component_dir);
                    bc.temperature = boundary.optional_temp;
                    bc.pressure = boundary.optional_pressure;
                    bc.boundary_faces.push_back({ boundary.i, boundary.j });
                    boundary.component_dir == 0 ? u_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open) : v_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open);
                }
                if ((boundary.component_dir == 0 && boundary.i != 0) || (boundary.component_dir == 1 && boundary.j != 0))
                {
                    auto& bc = open_boundary_conditions.emplace_back(1, boundary.component_dir);
                    bc.temperature = boundary.optional_temp;
                    bc.pressure = boundary.optional_pressure;
                    bc.boundary_faces.push_back({ boundary.i, boundary.j });
                    boundary.component_dir == 0 ? u_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open) : v_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open);
                }
            }
        }
    }

    std::vector<FixedInletBoundaryCondition> inlet_boundary_conditions;
    std::vector<FixedOutletBoundaryCondition> outlet_boundary_conditions;
    std::vector<OpenBoundaryCondition> open_boundary_conditions;
    FrictionBoundaryCondition friction_boundary_condition;

    //////////////
    // Settings //
    //////////////

    struct FluidProperties
    {
        double density{ 1.2 };
        double kinematic_viscosity{ 0.001 };
        double conductivity{ 0.1 };
        double cp{ 1000.0 };
    };
    FluidProperties fluid;

    int nx, ny;
    double dx, dy;
    double dt{ 10e-3 };
    double relaxation{ 0.5 };

    int innerPressureItterations{ 200 };
    int innerTemperatureItterations{ 20 };
    double residualLimit{ 10e-12 };

    double maxPressureResidual{ 0.0 }, maxDivergence{ 0.0 }, maxTemperatureResidual{ 0.0 };

    /////////////////
    // Solver Data //
    /////////////////

    // Cell Props

    PhysicalField p{ nx, ny };
    PhysicalField p_correction{ nx, ny };
    PhysicalField p_correction_old{ nx, ny };

    PhysicalField divergence{ nx, ny };
    PhysicalField p_scr{ nx, ny };
    PhysicalField p_coef{ nx, ny };

    PhysicalField t{ nx, ny };
    PhysicalField t_scr{ nx, ny };
    PhysicalField t_coeff{ nx, ny };

    PhysicalField residual{ nx, ny };
    PhysicalField t_previous{ nx, ny };

    // Face Props (staggered grid)

    PhysicalField u{ nx + 1, ny }, v{ nx, ny + 1 };

    PhysicalField u_scr{ nx + 1, ny };
    PhysicalField v_scr{ nx, ny + 1 };

    PhysicalField u_coeff{ nx + 1, ny };
    PhysicalField v_coeff{ nx, ny + 1 };

    FairMutex cell_data_mutex;

    PhysicalField cell_u{ nx, ny };
    PhysicalField cell_v{ nx, ny };

    // Face Flags (staggered grid): 1 - blocked, 2 - friction

    FlagField u_face_flags{ nx + 1, ny };
    FlagField v_face_flags{ nx, ny + 1 };
};