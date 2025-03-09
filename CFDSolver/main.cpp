#include "Solver.hpp"
#include "GUI/Interface.hpp"

int main()
{
    Interface interface;

    while (!WindowShouldClose())
    {
        // Define Model

        interface.RenderModel();

        Solver s(interface);

        // Solve for max n itterations

        auto solve_thread = std::jthread(&Solver::solve, &s, 1000000);

        interface.RenderResults();
        if (s.is_solving) // If it is still solving stop and post process
        {
            s.is_solving = false;
            solve_thread.join();
            interface.RenderResults();
        }
    }
	
	return 0;
}