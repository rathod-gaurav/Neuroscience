#pragma once

#include "OutputWriter.hpp"

template <unsigned int Nsd, unsigned int BfOrder>
class HodgkinHuxley{
    public:
        HodgkinHuxley(
            const double x_ll, const double x_ul, const unsigned int n_elements,
            const unsigned int NT, const double dt,
            const unsigned int quadOrder,
            const double Vrest, const double m_rest, const double h_rest, const double n_rest,
            const std::function<double(double)> &alpha_m, const std::function<double(double)> &beta_m,
            const std::function<double(double)> &alpha_h, const std::function<double(double)> &beta_h,
            const std::function<double(double)> &alpha_n, const std::function<double(double)> &beta_n,
            const double Cm, const double g_Na_max, const double g_K_max, const double g_L_max,
            const double V_Na, const double V_K, const double V_L, const double D,
            OutputWriter<Nsd,BfOrder>& output_writer
        );
        void run();
    
    private:
        const double x_ll_, x_ul_;
        const unsigned int n_elements_;
        const unsigned int NT_;
        const double dt_;
        const unsigned int quadOrder_;
        const double Vrest_, m_rest_, h_rest_, n_rest_;
        const std::function<double(double)> alpha_m_, beta_m_, alpha_h_, beta_h_, alpha_n_, beta_n_;
        const double Cm_, g_Na_max_, g_K_max_, g_L_max_;
        const double V_Na_, V_K_, V_L_, D_;
        OutputWriter<Nsd,BfOrder>& output_writer_;

        void make_grid();
        void setup_system(); 
        void apply_stimulus_current();
        void apply_stimulus_voltage();
        void compute_element(const typename DoFHandler<Nsd>::active_cell_iterator &elem, FEValues<Nsd>& fe_values, FullMatrix<double>& Mlocal, FullMatrix<double>& Klocal);
        void assemble_system();
        void compute_element_F_MG(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, FullMatrix<double>& MGlocal, Vector<double>& Flocal, std::vector<types::global_dof_index>& local_dof_indices);
        void assemble_system_F_MG();
        void solve();

        Triangulation<Nsd> triangulation;
        const FE_Q<Nsd> fe;
        DoFHandler<Nsd> dof_handler;
        unsigned int Nt; //total number of degrees of freedom
        std::map<types::global_dof_index, Point<Nsd>> dof_locations_map; //map of global dof :: point locations

        SparsityPattern sparsity_pattern;
        SparseMatrix<double> Mglobal; //global mass matrix 
        SparseMatrix<double> Kglobal; //global stiffness matrix
        SparseMatrix<double> MGglobal; //global ion-weighted mass matrix
        Vector<double> Fglobal; //global forcing vector + ion potential vector
        Vector<double> V, V_np1;
        Vector<double> m, m_np1;
        Vector<double> h, h_np1;
        Vector<double> n, n_np1;
        SparseMatrix<double> LHS;
        Vector<double> RHS;


        //Stimulus location and time
        const double x_stimulus_ll = 0.25; //location of stimulus
        const double x_stimulus_ul = 0.75; //location of stimulus
        const double t_stimulus_start = 1.0; //start time of stimulus current
        const double t_stimulus_end = 2.0; //end time of stimulus current //ms

        //Stimulus current
        Vector<double> Iext; //external stimulus current vector
        const double Iext_val = 20.0; //magnitude of stimulus current // micro-Amp/cm2
        
        //Stimulus voltage
        const double V_stimulus_val = -30.0; //magnitude of voltage perturbation // mV
};

#include "HodgkinHuxley.tpp"

#include "MakeGrid.tpp"
#include "SetupSystem.tpp"
#include "ApplyStimulusCurrent.tpp"
#include "ApplyStimulusVoltage.tpp"

#include "Solver.tpp"