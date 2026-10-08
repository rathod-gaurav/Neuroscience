#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void HodgkinHuxley<Nsd,BfOrder>::apply_stimulus_current(){
    for (unsigned int i = 0; i < Nt; ++i){
        if (x_stimulus_ll <= dof_locations_map[i][0] && dof_locations_map[i][0] <= x_stimulus_ul){
            Iext(i) = Iext_val;
        }
        else{
            Iext(i) = 0.0;
        }
    }
}