#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
HodgkinHuxley<Nsd,BfOrder>::HodgkinHuxley(
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
) : 
    x_ll_(x_ll), x_ul_(x_ul), n_elements_(n_elements),
    NT_(NT), dt_(dt), 
    quadOrder_(quadOrder),
    Vrest_(Vrest), m_rest_(m_rest), h_rest_(h_rest), n_rest_(n_rest),
    alpha_m_(alpha_m), beta_m_(beta_m), 
    alpha_h_(alpha_h), beta_h_(beta_h),
    alpha_n_(alpha_n), beta_n_(beta_n),
    Cm_(Cm), g_Na_max_(g_Na_max), g_K_max_(g_K_max), g_L_max_(g_L_max),
    V_Na_(V_Na), V_K_(V_K), V_L_(V_L), D_(D),
    output_writer_(output_writer),
    fe(BfOrder), 
    dof_handler(triangulation)
{}

template <unsigned int Nsd, unsigned int BfOrder>
void HodgkinHuxley<Nsd,BfOrder>::run(){
    make_grid();
    setup_system();
    solve();
}