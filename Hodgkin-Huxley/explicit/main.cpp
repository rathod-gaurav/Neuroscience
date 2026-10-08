// Hodgkin-Huxley fully explicit solver using deal.II library

#include <deal.II/grid/tria.h> //triangulation
#include <deal.II/dofs/dof_handler.h> //enumeration of degrees of freedom
#include <deal.II/grid/grid_generator.h> //grid generation

#include <deal.II/fe/fe_q.h> //Lagrange finite elements
#include <deal.II/dofs/dof_tools.h> //DoFHandler tools

#include <deal.II/fe/fe_values.h> //used to assemble matrix using quadrature on each cell
#include <deal.II/base/quadrature_lib.h> //quadrature rules

//need thiese three for treatment of boundary values
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

//for periodic boundary conditions
#include <deal.II/grid/grid_tools.h>

//linear algebra
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/base/table.h>
// #include <deal.II/lac/solver_cg.h>
// #include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/precondition.h>

//output
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <filesystem>

using namespace dealii;

#include <HodgkinHuxley.hpp>

int main(){
    constexpr unsigned int Nsd = 1; //2 for 2D, 3 for 3D
    constexpr unsigned int BfOrder = 1; //1 for linear, 2 for quadratic

    unsigned int quadOrder = 3; //quadrature order for numerical integration
    
    //Problem parameters
    const double x_ll = 0.0; //lower left corner of the domain
    const double x_ul = 10.0; //upper right corner of the domain
    const unsigned int n_elements = 1000; //number of global refinements in the mesh grid

    const unsigned int NT = 1000; //number of time steps
    const double dt = 0.1; //time step size

    const double Cm = 1.0; //membrane capacitance, uF/cm^2
    const double g_Na_max = 120.0; //maximum sodium conductance, mS/cm^2
    const double g_K_max = 36.0; //maximum potassium conductance, mS/cm^2
    const double g_L_max = 0.3; //leakage conductance, mS/cm^2

    const double V_Na = 50.0; //sodium reversal potential, mV
    const double V_K = -77.0; //potassium reversal potential, mV
    const double V_L = -54.4; //leakage reversal potential, mV

    const double a = 0.0238; //squid axon radius, cm
    const double R = 35.4; //squid axon resistance, ohm*cm
    const double D = (1e3*a)/(2.0*R); //diffusion coefficient, mS
    // std::cout << "D = " << D << std::endl;

    auto alpha_m = [](const double V) { return (0.1 * (V + 40.0)) / (1.0 - std::exp(-0.1 * (V + 40.0))); };
    auto beta_m = [](const double V) { return 4.0 * std::exp(-(V + 65.0)/18.0); };
    auto alpha_h = [](const double V) { return 0.07 * std::exp(-(V + 65.0)/20.0); };
    auto beta_h = [](const double V) { return 1.0 / (1.0 + std::exp(-0.1 * (V + 35.0))); };
    auto alpha_n = [](const double V) { return (0.01 * (V + 55.0)) / (1.0 - std::exp(-0.1 * (V + 55.0))); };
    auto beta_n = [](const double V) { return 0.125 * std::exp(-(V + 65.0)/80.0); };

    const double Vrest = -65.0; //resting potential, mV

    double m_inf = alpha_m(Vrest) / (alpha_m(Vrest) + beta_m(Vrest));
    double h_inf = alpha_h(Vrest) / (alpha_h(Vrest) + beta_h(Vrest));
    double n_inf = alpha_n(Vrest) / (alpha_n(Vrest) + beta_n(Vrest));

    const double m_rest = m_inf; //resting value of m //should be 0.053
    const double h_rest = h_inf; //resting value of h //should be 0.596
    const double n_rest = n_inf; //resting value of n //should be 0.318

    // std::cout << "Resting values: m_rest = " << m_rest << ", h_rest = " << h_rest << ", n_rest = " << n_rest << std::endl;

    OutputWriter<Nsd,BfOrder> output_writer("output");
    HodgkinHuxley<Nsd,BfOrder> problem(
        x_ll, x_ul, n_elements,
        NT, dt,
        quadOrder,
        Vrest, m_rest, h_rest, n_rest,
        alpha_m, beta_m,
        alpha_h, beta_h,
        alpha_n, beta_n,
        Cm, g_Na_max, g_K_max, g_L_max,
        V_Na, V_K, V_L, D,
        output_writer
    );
    
    problem.run();
}