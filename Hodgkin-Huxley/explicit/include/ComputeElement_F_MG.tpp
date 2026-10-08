#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void HodgkinHuxley<Nsd,BfOrder>::compute_element_F_MG(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, FullMatrix<double>& MGlocal, Vector<double>& Flocal, std::vector<types::global_dof_index>& local_dof_indices){
    fe_values.reinit(elem);

    for(const unsigned int q_index : fe_values.quadrature_point_indices()){
        double Iext_local = 0.0;
        double m_np1_local = 0.0;
        double h_np1_local = 0.0;
        double n_np1_local = 0.0;
        for(const unsigned int A : fe_values.dof_indices()){
            Iext_local += fe_values.shape_value(A,q_index)*Iext(local_dof_indices[A]);
            m_np1_local += fe_values.shape_value(A,q_index)*m_np1(local_dof_indices[A]);
            h_np1_local += fe_values.shape_value(A,q_index)*h_np1(local_dof_indices[A]);
            n_np1_local += fe_values.shape_value(A,q_index)*n_np1(local_dof_indices[A]);
        }

        // G_np1_local = g_Na_max_*std::pow(m_np1_local,3)*h_np1_local + g_K_max_*std::pow(n_np1_local,4) + g_L_max_;
        double g_Na_q = g_Na_max_*std::pow(m_np1_local,3)*h_np1_local;
        double g_K_q = g_K_max_*std::pow(n_np1_local,4);
        double g_L_q = g_L_max_;
        double G_np1_local = g_Na_q + g_K_q + g_L_q;
        double E_np1_local = g_Na_q*V_Na_ + g_K_q*V_K_ + g_L_q*V_L_;

        for(const unsigned int A : fe_values.dof_indices()){
            for(const unsigned int B : fe_values.dof_indices()){
                MGlocal(A,B) += (G_np1_local*fe_values.shape_value(A,q_index)*fe_values.shape_value(B,q_index)*fe_values.JxW(q_index));
            }

            Flocal(A) += ((Iext_local + E_np1_local)*fe_values.shape_value(A,q_index)*fe_values.JxW(q_index));
        }
    }
}