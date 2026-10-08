#pragma once

#include "AssembleSystem.tpp"
#include "AssembleSystem_F_MG.tpp"

template <unsigned int Nsd, unsigned int BfOrder>
void HodgkinHuxley<Nsd,BfOrder>::solve(){
    std::cout << "Initiating solver..." << std::endl;

    std::cout << "-----------------------------" << std::endl;
    Mglobal = 0.0;
    Kglobal = 0.0;
    assemble_system();

    output_writer_.write_vtu(dof_handler, V, "V", m, "m", h, "h", n, "n", 0); //write initial condition to output file
    // Iext = 0.0;

    double t = dt_;
    for(unsigned int timestep = 1; timestep < NT_; ++timestep){

        
        // if(t >= t_stimulus_start && t <= t_stimulus_end){
        //     apply_stimulus_current();
        // } 
        // else{
        //     Iext = 0.0;
        // }    
        
        //constant voltage stimulus
        // apply_stimulus_voltage();    
        
        //constant current stimulus
        if(t >= t_stimulus_start){
            apply_stimulus_current();
        }
        
        

        //Solve ODEs first to get updated gating variables - finite difference method
        for(unsigned int i = 0; i < Nt; ++i){
            double V_i = V(i);
            double m_i = m(i);
            double h_i = h(i);
            double n_i = n(i);

            // std::cout << "V_i = " << V_i << ", m_i = " << m_i << ", h_i = " << h_i << ", n_i = " << n_i << std::endl;

            m_np1(i) = m_i + dt_*(alpha_m_(V_i) * (1.0 - m_i) - beta_m_(V_i) * m_i);
            h_np1(i) = h_i + dt_*(alpha_h_(V_i) * (1.0 - h_i) - beta_h_(V_i) * h_i);
            n_np1(i) = n_i + dt_*(alpha_n_(V_i) * (1.0 - n_i) - beta_n_(V_i) * n_i);

            // std::cout << "m_np1(i) = " << m_np1(i) << ", h_np1(i) = " << h_np1(i) << ", n_np1(i) = " << n_np1(i) << std::endl;

        }

        //Now solve PDE to get updated voltage
        Fglobal = 0.0;
        MGglobal = 0.0;
        assemble_system_F_MG();
        
        // LHS = Cm*Mglobal + dt*D*Kglobal + dt*MGglobal;
        LHS = 0.0;
        LHS.add(Cm_,Mglobal);
        LHS.add(dt_*D_,Kglobal);
        LHS.add(dt_,MGglobal);

        //RHS = dt*Fglobal + Cm*Mglobal*V;
        RHS = 0.0;
        Mglobal.vmult(RHS,V);
        RHS *= Cm_;
        RHS.add(dt_,Fglobal);

        //Solve for V_np1
        SolverControl control(1000, 1e-10*RHS.l2_norm());
        SolverCG<Vector<double>> cgsolver(control);
        cgsolver.solve(LHS, V_np1, RHS, PreconditionIdentity());

        // std::cout << "----------" << std::endl;
        // std::cout << "V_Na = " << g_Na_max_*m_np1(0)*m_np1(0)*m_np1(0)*h_np1(0)*(V_np1(0) - V_Na_) << std::endl;
        // std::cout << "V_K = " << g_K_max_*n_np1(0)*n_np1(0)*n_np1(0)*n_np1(0)*(V_np1(0) - V_K_) << std::endl;
        // std::cout << "V_L = " << g_L_max_*(V_np1(0) - V_L_) << std::endl;
        // std::cout << "----------" << std::endl;
        
        V = V_np1;
        m = m_np1;
        h = h_np1;
        n = n_np1;
        t += dt_;
        std::cout << "Solve completed for timestep: " << timestep << std::endl;\
        
        output_writer_.write_vtu(dof_handler, V, "V", m, "m", h, "h", n, "n", timestep);
    }

    std::cout << "Solve completed. Congratulations!!" << std::endl;

    output_writer_.write_pvd();
}