#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void HodgkinHuxley<Nsd,BfOrder>::setup_system(){
    dof_handler.distribute_dofs(fe);
    Nt = dof_handler.n_dofs();

    std::cout << "Number of degrees of freedom: " << Nt << std::endl;

    DynamicSparsityPattern dsp(Nt, Nt);
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    Mglobal.reinit(sparsity_pattern);
    Kglobal.reinit(sparsity_pattern);
    MGglobal.reinit(sparsity_pattern);
    Fglobal.reinit(Nt);

    V.reinit(Nt);
    V_np1.reinit(Nt);
    m.reinit(Nt);
    m_np1.reinit(Nt);
    h.reinit(Nt);
    h_np1.reinit(Nt);
    n.reinit(Nt);
    n_np1.reinit(Nt);

    LHS.reinit(sparsity_pattern);
    RHS.reinit(Nt);

    //initial conditions
    V = Vrest_;
    m = m_rest_;
    h = h_rest_;
    n = n_rest_;
    
    dof_locations_map = DoFTools::map_dofs_to_support_points(MappingQ1<Nsd>(), dof_handler);

    //stimulus current vector
    Iext.reinit(Nt);
}