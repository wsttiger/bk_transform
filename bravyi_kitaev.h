/*******************************************************************************
 * Copyright (c) 2024 NVIDIA Corporation & Affiliates.                         *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/
#pragma once

// #include "cudaq/solvers/operators/molecule/fermion_compiler.h"

#include "spin_op.h"
#include "tensor.h"

namespace cudaq::solvers {

/// @brief Helper function used by the Bravyi-Kitaev transformation.
cudaq::spin_op seeley_richard_love(std::size_t i, std::size_t j, std::complex<double> coef, int n_qubits);

/// @brief Map fermionic operators to spin operators via the
/// Bravyi-Kitaev transformation.
// class bravyi_kitaev : public fermion_compiler {
class bravyi_kitaev {
public:
//   cudaq::spin_op generate(const double constant, const cudaqx::tensor<> &hpq,
//                           const cudaqx::tensor<> &hpqrs) override;
  cudaq::spin_op generate(const double constant, const cudaqx::tensor<> &hpq,
                          const cudaqx::tensor<> &hpqrs);

//  CUDAQ_EXTENSION_CREATOR_FUNCTION(fermion_compiler, jordan_wigner)
};
// CUDAQ_REGISTER_TYPE(jordan_wigner)
} // namespace cudaq::solvers
