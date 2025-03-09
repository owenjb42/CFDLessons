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
        /* TO-DO */
    }

    void ComputeConvectiveTerms()
    {
        /* TO-DO */
    }

    void ComputeDiffusiveTerms()
    {
        /* TO-DO */
    }

    void ComputeConvectiveFluxesAndDiffusiveTermsForTemperature()
    {
        /* TO-DO */
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
        /* TO-DO */
    }

    void ApplyPressureCorrectionToVelocity()
    {
        /* TO-DO */
    }

    void SolveVelocitiesForMomentumEquation()
    {
        /* TO-DO */
    }

    void SolveTemperatures()
    {
        /* TO-DO */
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
    }

    void ApplyPressureBoundaryConditions()
    {
        for (auto& boundary : open_boundary_conditions)
            boundary.ApplyForPressure(*this);
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