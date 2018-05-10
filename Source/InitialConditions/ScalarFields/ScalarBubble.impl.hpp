/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARBUBBLE_HPP_)
#error "This file should only be included through ScalarBubble.hpp"
#endif

#ifndef SCALARBUBBLE_IMPL_HPP_
#define SCALARBUBBLE_IMPL_HPP_

inline ScalarBubble::ScalarBubble(params_t a_params, double a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

// Compute the value of the initial vars on the grid
template <class data_t>
void ScalarBubble::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ScalarField<>>::Vars<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.centerSF);

    data_t xx = coords.x;
    data_t yy = coords.y;
    data_t zz = coords.z;
    data_t rr = coords.get_radius();

    data_t Input1;
    data_t Input2; 

    double m_spacing = m_params.spacing;
    int indxL = static_cast<int>(floor(rr/m_spacing));
    int indxH = static_cast<int>(ceil(rr/m_spacing));

    Input1 =  m_params.inputValues1[indxL]
                + (rr/m_spacing - indxL)*(m_params.inputValues1[indxH] - m_params.inputValues1[indxL]);
    Input2 =  m_params.inputValues2[indxL]
                + (rr/m_spacing - indxL)*(m_params.inputValues2[indxH] - m_params.inputValues2[indxL]);
    // set the field vars
    vars.phi = 0;  
    vars.Pi = Input1;

    // start with unit lapse and flat metric (must be relaxed for chi)
    vars.lapse = 1;
    vars.chi = 1;

    // conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

// Compute the value of phi at the current point
template <class data_t>
data_t ScalarBubble::compute_phi(Coordinates<data_t> coords) const
{
    data_t rr = coords.get_radius();
    data_t rr2 = rr * rr;
    data_t out_phi = m_params.amplitudeSF * rr2 *
                     exp(-pow(rr - m_params.r_zero / m_params.widthSF, 2.0));

    return out_phi;
}

#endif /* SCALARBUBBLE_IMPL_HPP_ */
