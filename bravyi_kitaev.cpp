#include <set>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cassert>

#include "bravyi_kitaev.h"

namespace cudaq::solvers {

void print_set(const std::set<std::size_t>& s) {
    if (s.empty()) {
        std::cout << "[]\n";
        return;
    }
    std::cout << "[";
    std::copy(s.begin(), s.end(), 
        std::ostream_iterator<std::size_t>(std::cout, ", "));
    std::cout << "\b\b]\n"; // backspace over last comma and space
}

template<typename T>
std::set<T> set_difference(const std::set<T>& set1, const std::set<T>& set2) {
    std::set<T> result;
    std::set_difference(set1.begin(),
                        set1.end(),
                        set2.begin(),
                        set2.end(),
                        std::inserter(result, result.begin()));
    return result;
}

template<typename T>
std::set<T> set_union(const std::set<T>& set1, const std::set<T>& set2) {
    std::set<T> result;
    std::set_union(set1.begin(),
                   set1.end(),
                   set2.begin(),
                   set2.end(),
                   std::inserter(result, result.begin()));
    return result;
}

template<typename T>
std::set<T> set_intersection(const std::set<T>& set1, const std::set<T>& set2) {
    std::set<T> result;
    std::set_intersection(set1.begin(),
                          set1.end(),
                          set2.begin(),
                          set2.end(),
                          std::inserter(result, result.begin()));
    return result;
}

template<typename T>
std::set<T> set_symmetric_difference(const std::set<T>& set1, const std::set<T>& set2) {
    std::set<T> result;
    std::set_symmetric_difference(set1.begin(),
                                  set1.end(), 
                                  set2.begin(), 
                                  set2.end(),
                                  std::inserter(result, result.begin()));
    return result;
}

std::set<std::size_t> occupation_set(std::size_t index) {
    std::set<std::size_t> indices;
    index += 1;
    indices.insert(index - 1);
    
    std::size_t parent = index & (index - 1);
    index -= 1;
    
    while (index != parent) {
        indices.insert(index - 1);
        index &= index - 1;
    }
    
    return indices;
}

std::set<std::size_t> parity_set(std::size_t index) {
    std::set<std::size_t> indices;
    
    while (index > 0) {
        indices.insert(index - 1);
        index &= index - 1;
    }
    
    return indices;
}

std::set<std::size_t> update_set(std::size_t index, std::size_t n_qubits) {
    std::set<std::size_t> indices;
    
    index += 1;
    index += index & -index;
    
    while (index <= n_qubits) {
        indices.insert(index - 1);
        index += index & -index;
    }
    
    return indices;
}

std::set<std::size_t> remainder_set(int index) {
    return set_difference(parity_set(index), occupation_set(index));
}

std::set<std::size_t> F_ij_set(std::size_t i, std::size_t j) {
    return set_symmetric_difference(occupation_set(i), occupation_set(j));
}

std::set<std::size_t> P0_ij_set(std::size_t i, std::size_t j) {
    return set_symmetric_difference(parity_set(i), parity_set(j));
}

std::set<std::size_t> P1_ij_set(std::size_t i, std::size_t j) {
    return set_symmetric_difference(parity_set(i), remainder_set(j));
}

std::set<std::size_t> P2_ij_set(std::size_t i, std::size_t j) {
    return set_symmetric_difference(remainder_set(i), parity_set(j));
}

std::set<std::size_t> P3_ij_set(std::size_t i, std::size_t j) {
    return set_symmetric_difference(remainder_set(i), remainder_set(j));
}

std::set<std::size_t> U_ij_set(std::size_t i, std::size_t j, std::size_t n_qubits) {
    return set_symmetric_difference(update_set(i, n_qubits), update_set(j, n_qubits));
}

std::set<std::size_t> alpha_set(std::size_t i, std::size_t j, std::size_t n_qubits) {
    return set_intersection(update_set(i, n_qubits), parity_set(j));
}

std::set<std::size_t> U_diff_a_set(std::size_t i, std::size_t j, std::size_t n_qubits) {
    return set_difference(U_ij_set(i, j, n_qubits), alpha_set(i, j, n_qubits));
}

std::set<std::size_t> P0_ij_diff_a_set(std::size_t i, std::size_t j, std::size_t n_qubits) {
    return set_symmetric_difference(P0_ij_set(i, j), alpha_set(i, j, n_qubits));
}

cudaq::spin_op seeley_richard_love(std::size_t i, std::size_t j, std::complex<double> coef, int n_qubits) {

    using double_complex = std::complex<double>;
    const double_complex imag_i = double_complex(0.0, 1.0);

    coef *= 0.25;

    cudaq::spin_op seeley_richard_love_result = 0.0 * cudaq::spin::i(0);

    // Case 0
    if (i == j) {
        std::cout << "Case 0\n";
        if (not occupation_set(i).empty()) {
            cudaq::spin_op ops;
            for (int index : occupation_set(i)) {
                ops *= cudaq::spin::z(index);
            }
            seeley_richard_love_result += -coef * 2.0 * ops;
        }
        seeley_richard_love_result +=  coef * 2.0 * cudaq::spin::i(0);
    }

    // Case 1
    else if (i % 2 == 0 and j % 2 == 0) {
        std::cout << "Case 1\n";
        cudaq::spin_op x_pad;
        for (int index : U_diff_a_set(i, j, n_qubits)) {
            x_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op y_pad;
        for (int index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }
        cudaq::spin_op z_pad;
        for (int index : P0_ij_diff_a_set(i, j, n_qubits)) {
            z_pad *= cudaq::spin::z(index);
        }
        cudaq::spin_op left_pad = x_pad * y_pad * z_pad;

        cudaq::spin_op op1 = left_pad * cudaq::spin::y(j) * cudaq::spin::x(i);
        cudaq::spin_op op2 = left_pad * cudaq::spin::x(j) * cudaq::spin::y(i);
        cudaq::spin_op op3 = left_pad * cudaq::spin::x(j) * cudaq::spin::x(i);
        cudaq::spin_op op4 = left_pad * cudaq::spin::y(j) * cudaq::spin::y(i);

        if (i < j) {
            seeley_richard_love_result +=           coef * op1;
            seeley_richard_love_result +=          -coef * op2;
            seeley_richard_love_result += -imag_i * coef * op3;
            seeley_richard_love_result += -imag_i * coef * op4;
        } else {  
            seeley_richard_love_result += -imag_i * coef * op1;
            seeley_richard_love_result +=  imag_i * coef * op2;
            seeley_richard_love_result +=          -coef * op3;
            seeley_richard_love_result +=          -coef * op4;
        }
    }

    // Case 2
    else if (i % 2 == 1 and j % 2 == 0 and not parity_set(j).contains(i)) {
        std::cout << "Case 2\n";
        cudaq::spin_op x_pad;
        for (int index : U_diff_a_set(i, j, n_qubits)) {
            x_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op y_pad;
        for (int index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }
        cudaq::spin_op left_pad = x_pad * y_pad;

        cudaq::spin_op right_pad_1;
        auto P0_minus_alpha = set_difference(P0_ij_set(i, j), alpha_set(i, j, n_qubits));
        for (int index : P0_minus_alpha) {
            right_pad_1 *= cudaq::spin::z(index);
        }

        cudaq::spin_op right_pad_2;
        auto P2_minus_alpha = set_difference(P2_ij_set(i, j), alpha_set(i, j, n_qubits));
        for (int index : P2_minus_alpha) {
            right_pad_2 *= cudaq::spin::z(index);
        }

        double_complex c0, c1, c2, c3;
        if (i < j) {
            c0 =           coef;
            c1 = -imag_i * coef;
            c2 =          -coef;
            c3 = -imag_i * coef;
        } else {
            c0 = -imag_i * coef;
            c1 =          -coef;
            c2 =  imag_i * coef;
            c3 =          -coef;
        }

        seeley_richard_love_result += c0 * left_pad * cudaq::spin::y(j) * cudaq::spin::x(i) * right_pad_1;
        seeley_richard_love_result += c1 * left_pad * cudaq::spin::x(j) * cudaq::spin::x(i) * right_pad_1;
        seeley_richard_love_result += c2 * left_pad * cudaq::spin::x(j) * cudaq::spin::y(i) * right_pad_2;
        seeley_richard_love_result += c3 * left_pad * cudaq::spin::y(j) * cudaq::spin::y(i) * right_pad_2;
    }

    // Case 3
    else if (i % 2 == 1 and j % 2 == 0 and parity_set(j).contains(i)) {
        std::cout << "Case 3\n";
        cudaq::spin_op left_pad;
        for (auto index : U_ij_set(i, j, n_qubits)) {
            left_pad *= cudaq::spin::x(index);
        }

        cudaq::spin_op right_pad_1;
        for (auto index : set_difference(P0_ij_set(i, j), {i})) {
            right_pad_1 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_2;
        for (auto index : set_difference(P2_ij_set(i, j), {i})) {
            right_pad_2 *= cudaq::spin::z(index);
        }

        seeley_richard_love_result +=           coef * left_pad * cudaq::spin::y(j) * cudaq::spin::y(i) * right_pad_1;
        seeley_richard_love_result += -imag_i * coef * left_pad * cudaq::spin::x(j) * cudaq::spin::y(i) * right_pad_1;
        seeley_richard_love_result +=           coef * left_pad * cudaq::spin::x(j) * cudaq::spin::x(i) * right_pad_2;
        seeley_richard_love_result +=  imag_i * coef * left_pad * cudaq::spin::y(j) * cudaq::spin::x(i) * right_pad_2;
    }

    // Case 4
    else if (i % 2 == 0 and j % 2 == 1 and not parity_set(j).contains(i) and not update_set(i, n_qubits).contains(i)) {
        std::cout << "Case 4\n";
        cudaq::spin_op x_pad;
        for (auto index : U_diff_a_set(i, j, n_qubits)) {
            x_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op y_pad;
        for (auto index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }
        auto left_pad = x_pad * y_pad;

        cudaq::spin_op right_pad_1;
        for (auto index : set_difference(P0_ij_set(i, j), alpha_set(i, j, n_qubits))) {
            right_pad_1 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_2;
        for (auto index : set_difference(P1_ij_set(i, j), alpha_set(i, j, n_qubits))) {
            right_pad_2 *= cudaq::spin::z(index);
        }
     
        double_complex c0, c1, c2, c3;
        if (i < j) {
            c0 =          -coef;
            c1 = -imag_i * coef;
            c2 =           coef;
            c3 = -imag_i * coef;
        } else {
            c0 =  imag_i * coef;
            c1 =          -coef;
            c2 = -imag_i * coef;
            c3 =          -coef;
        }
        seeley_richard_love_result += c0 * left_pad * cudaq::spin::x(j) * cudaq::spin::y(i) * right_pad_1;
        seeley_richard_love_result += c1 * left_pad * cudaq::spin::x(j) * cudaq::spin::x(i) * right_pad_1;
        seeley_richard_love_result += c2 * left_pad * cudaq::spin::y(j) * cudaq::spin::x(i) * right_pad_2;
        seeley_richard_love_result += c3 * left_pad * cudaq::spin::y(j) * cudaq::spin::y(i) * right_pad_2;
    }

    // Case 5
    else if (i % 2 == 0 and j % 2 == 1 and not parity_set(j).contains(i) and update_set(i, n_qubits).contains(j)) {
        std::cout << "Case 5\n";
        cudaq::spin_op left_pad_1;
        auto x_range_1 = set_difference(U_ij_set(i, j, n_qubits), {j});
        for (auto index : x_range_1) {
            left_pad_1 *= cudaq::spin::x(index);
        }
        cudaq::spin_op left_pad_2;
        auto x_range_2 = set_difference(x_range_1, alpha_set(i, j, n_qubits));
        for (auto index : x_range_2) {
            left_pad_2 *= cudaq::spin::x(index);
        }

        cudaq::spin_op y_pad;
        for (auto index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }
        cudaq::spin_op z_pad;
        for (auto index : set_difference(P0_ij_set(i, j), alpha_set(i, j, n_qubits))) {
            z_pad *= cudaq::spin::z(index);
        }
        auto right_pad_1 = y_pad * z_pad;

        cudaq::spin_op right_pad_2;
        for (auto index : set_union(P1_ij_set(i, j), {j})) {
            right_pad_2 *= cudaq::spin::z(index);
        }
    }

    // Case 6
    else if (i % 2 == 0 and j % 2 == 1 and parity_set(j).contains(i) and update_set(i, n_qubits).contains(j)) {
        std::cout << "Case 6\n";
        cudaq::spin_op left_pad;
        for (auto index : set_difference(U_ij_set(i, j, n_qubits), {j})) {
            left_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op right_pad;
        for (auto index : set_union(P1_ij_set(i, j), {j})) {
            right_pad *= cudaq::spin::z(index);
        }

        seeley_richard_love_result +=           coef * left_pad * cudaq::spin::x(i);
        seeley_richard_love_result += -imag_i * coef * left_pad * cudaq::spin::y(i);
        seeley_richard_love_result +=  imag_i * coef * left_pad * cudaq::spin::y(i) * right_pad;
        seeley_richard_love_result +=          -coef * left_pad * cudaq::spin::x(i) * right_pad;
    }

    // Case 7
    else if (i % 2 == 1 and j % 2 == 1 and not parity_set(j).contains(i) and not update_set(i, n_qubits).contains(j)) {
        std::cout << "Case 7\n";
        cudaq::spin_op x_pad;
        for (auto index : U_diff_a_set(i, j, n_qubits)) {
            x_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op y_pad;
        for (auto index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }
        auto left_pad = x_pad * y_pad;

        cudaq::spin_op right_pad_1;
        for (auto index : set_difference(P0_ij_set(i, j), alpha_set(i, j, n_qubits))) {
            right_pad_1 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_2;
        for (auto index : set_difference(P1_ij_set(i, j), alpha_set(i, j, n_qubits))) {
            right_pad_2 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_3;
        for (auto index : set_difference(P2_ij_set(i, j), alpha_set(i, j, n_qubits))) {
            right_pad_3 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_4;
        for (auto index : set_difference(P3_ij_set(i, j), alpha_set(i, j, n_qubits))) {
            right_pad_4 *= cudaq::spin::z(index);
        }

        double_complex c0, c1, c2, c3;
        if (i < j) {
            c0 = -imag_i * coef;
            c1 =           coef;
            c2 =          -coef;
            c3 = -imag_i * coef;
        } else {
            c0 =          -coef;
            c1 = -imag_i * coef;
            c2 =  imag_i * coef;
            c3 =          -coef;
        }

        seeley_richard_love_result += c0 * left_pad * cudaq::spin::x(j) * cudaq::spin::x(i) * right_pad_1;
        seeley_richard_love_result += c1 * left_pad * cudaq::spin::y(j) * cudaq::spin::x(i) * right_pad_2;
        seeley_richard_love_result += c2 * left_pad * cudaq::spin::x(j) * cudaq::spin::y(i) * right_pad_3;
        seeley_richard_love_result += c3 * left_pad * cudaq::spin::y(j) * cudaq::spin::y(i) * right_pad_4;
    }

    // Case 8
    else if (i % 2 == 1 and j % 2 == 1 and parity_set(j).contains(i) and not update_set(i, n_qubits).contains(j)) {
        std::cout << "Case 8\n";
        cudaq::spin_op left_pad;
        for (auto index : U_ij_set(i, j, n_qubits)) {
            left_pad *= cudaq::spin::x(index);
        }

        cudaq::spin_op right_pad_1;
        for (auto index : set_difference(P0_ij_set(i, j), {i})) {
            right_pad_1 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_2;
        for (auto index : set_difference(P1_ij_set(i, j), {i})) {
            right_pad_2 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_3;
        for (auto index : set_difference(P2_ij_set(i, j), {i})) {
            right_pad_3 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_4;
        for (auto index : set_difference(P3_ij_set(i, j), {i})) {
            right_pad_4 *= cudaq::spin::z(index);
        }

        seeley_richard_love_result += -imag_i * coef * left_pad * cudaq::spin::x(j) * cudaq::spin::y(i) * right_pad_1;
        seeley_richard_love_result +=           coef * left_pad * cudaq::spin::y(j) * cudaq::spin::y(i) * right_pad_2;
        seeley_richard_love_result +=           coef * left_pad * cudaq::spin::x(j) * cudaq::spin::x(i) * right_pad_3;
        seeley_richard_love_result +=  imag_i * coef * left_pad * cudaq::spin::y(j) * cudaq::spin::x(i) * right_pad_4;
    }

    // Case 9
    else if (i % 2 == 1 and j % 2 == 1 and not parity_set(j).contains(i) and update_set(i, n_qubits).contains(j)) {
        std::cout << "Case 9\n";
        cudaq::spin_op left_pad_3;
        auto x_range_1 = set_difference(U_ij_set(i, j, n_qubits), {j}); 
        for (auto index : x_range_1) {
            left_pad_3 *= cudaq::spin::x(index);
        }

        cudaq::spin_op left_pad_1;
        auto x_range_2 = set_difference(x_range_1, alpha_set(i, j, n_qubits));
        for (auto index : x_range_2) {
            left_pad_1 *= cudaq::spin::x(index);
        }

        cudaq::spin_op left_pad_2;
        auto x_range_3 = set_difference(set_union(x_range_1, {i}), alpha_set(i, j, n_qubits));
        for (auto index : x_range_3) {
            left_pad_2 *= cudaq::spin::x(index);
        }

        auto z_range_1 = set_difference(P2_ij_set(i, j), alpha_set(i, j, n_qubits));
        cudaq::spin_op z_pad_1;
        for (auto index : z_range_1) {
            z_pad_1 *= cudaq::spin::z(index);
        }
        cudaq::spin_op y_pad;
        for (auto index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }
        auto right_pad_1 = z_pad_1 * y_pad;

        auto z_range_2 = set_difference(P0_ij_set(i, j), alpha_set(i, j, n_qubits));
        cudaq::spin_op z_pad_2;
        for (auto index : z_range_2) {
            z_pad_2 *= cudaq::spin::z(index);
        }
        auto right_pad_2 = z_pad_2 * y_pad;
        
        auto z_range_3 = set_union(P1_ij_set(i, j), {j});
        cudaq::spin_op right_pad_3;
        for (auto index : z_range_3) {
            right_pad_3 *= cudaq::spin::z(index);
        }

        auto z_range_4 = set_union(P3_ij_set(i, j), {j});
        cudaq::spin_op right_pad_4;
        for (auto index : z_range_4) {
            right_pad_4 *= cudaq::spin::z(index);
        }

        seeley_richard_love_result +=          -coef * left_pad_1 * cudaq::spin::y(i) * right_pad_1;
        seeley_richard_love_result += -imag_i * coef * left_pad_2 * right_pad_2;
        seeley_richard_love_result +=          -coef * left_pad_3 * cudaq::spin::x(i) * right_pad_3;
        seeley_richard_love_result +=  imag_i * coef * left_pad_3 * cudaq::spin::y(i) * right_pad_4;
    }

    // Case 10
    else if (i % 2 == 1 and j % 2 == 1 and parity_set(j).contains(i) and update_set(i, n_qubits).contains(j)) {
        std::cout << "Case 10\n";
        cudaq::spin_op left_pad;
        for (auto index : set_difference(U_ij_set(i, j, n_qubits), {j})) {
            left_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op right_pad_1;
        for (auto index : set_difference(P0_ij_set(i, j), {i})) {
            right_pad_1 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_2;
        for (auto index : set_difference(P2_ij_set(i, j), {i})) {
            right_pad_2 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_3;
        for (auto index : P1_ij_set(i, j)) {
            right_pad_3 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_4;
        for (auto index : P3_ij_set(i, j)) {
            right_pad_4 *= cudaq::spin::z(index);
        }

        seeley_richard_love_result += -imag_i * coef * left_pad * cudaq::spin::y(i) * right_pad_1;
        seeley_richard_love_result +=           coef * left_pad * cudaq::spin::x(i) * right_pad_2;
        seeley_richard_love_result +=          -coef * left_pad * cudaq::spin::z(j) * cudaq::spin::x(i) * right_pad_3;
        seeley_richard_love_result +=  imag_i * coef * left_pad * cudaq::spin::z(j) * cudaq::spin::y(i) * right_pad_4;
    }
    else {
        std::cout << "uh oh\n";
    }

    return seeley_richard_love_result;
}

cudaq::spin_op hermitian_one_body_product(std::size_t a, 
                                          std::size_t b, 
                                          std::size_t c, 
                                          std::size_t d, 
                                          std::complex<double> coef, 
                                          std::size_t nqubits) {

    cudaq::spin_op c_dag_c_ac = seeley_richard_love(a, c, coef, nqubits);
    cudaq::spin_op c_dag_c_bd = seeley_richard_love(b, d, 1, nqubits);
    c_dag_c_ac *= c_dag_c_bd;

    cudaq::spin_op hermitian_sum = c_dag_c_ac;

    cudaq::spin_op c_dag_c_ca = seeley_richard_love(c, a, std::conj(coef), nqubits);
    cudaq::spin_op c_dag_c_db = seeley_richard_love(d, b, 1, nqubits);
    c_dag_c_ca *= c_dag_c_db;

    hermitian_sum += c_dag_c_ca;
    return hermitian_sum;
}

std::complex<double> two_body_coef(const cudaqx::tensor<> &hpqrs, 
                                   std::size_t p,
                                   std::size_t q,
                                   std::size_t r,
                                   std::size_t s) {
    return hpqrs.at({p, q, r, s}) - hpqrs.at({q, p, r, s}) - hpqrs.at({p, q, s, r}) + hpqrs.at({q, p, s, r});
}

cudaq::spin_op bravyi_kitaev::generate(const double constant,
                                       const cudaqx::tensor<> &hpq,
                                       const cudaqx::tensor<> &hpqrs) {
    assert(hpq.rank() == 2 && "hpq must be a rank-2 tensor");
    assert(hpqrs.rank() == 4 && "hpqrs must be a rank-4 tensor");
    auto nqubits = hpq.shape()[0];
    cudaq::spin_op bk_hamiltonian = 0.0 * cudaq::spin::i(0);
    double constant_term = constant;
    for (std::size_t p = 0; p < nqubits; p++) {
        if (std::abs(hpq.at({p,p})) > 0.0) {
            bk_hamiltonian += seeley_richard_love(p, p, hpq.at({p,p}), nqubits); 
        }
        for (std::size_t q = 0; q < p; q++) {
            if (std::abs(hpq.at({p,q})) > 0.0) {
                bk_hamiltonian += seeley_richard_love(p, q,           hpq.at({p,q} ), nqubits); 
                bk_hamiltonian += seeley_richard_love(q, p, std::conj(hpq.at({p,q})), nqubits); 
            }
            auto coef = 0.25 * two_body_coef(hpqrs, p, q, q, p);
            if (std::abs(coef) > 0.0) {
                if (not occupation_set(p).empty()) {
                    cudaq::spin_op zs;
                    for (auto index : occupation_set(p)) {
                        zs *= cudaq::spin::z(index);
                    }
                    bk_hamiltonian += -coef * zs;
                }
                if (not occupation_set(q).empty()) {
                    cudaq::spin_op zs2;
                    for (auto index : occupation_set(q)) {
                        zs2 *= cudaq::spin::z(index);
                    }
                    bk_hamiltonian += -coef * zs2;
                }
                if (not F_ij_set(p, q).empty()) {
                    cudaq::spin_op zs3;
                    for (auto index : F_ij_set(p, q)) {
                        zs3 *= cudaq::spin::z(index);
                    }
                    bk_hamiltonian +=  coef * zs3;
                }
                // TODO: question about this:
                // constant_term is a double but += will cast it to a std::complex<>
                // 
                // constant_term += coef;
                constant_term += std::real(coef);
            }
        }
    }

    for (std::size_t p = 0; p < nqubits; p++) {
        for (std::size_t q = 0; q < nqubits; q++) {
            for (std::size_t r = 0; r < q; r++) {
                if ((p != q) and (p != r)) {
                    auto coef = two_body_coef(hpqrs, p, q, r, p);
                    if (std::abs(coef) > 0.0) {
                        cudaq::spin_op excitation = seeley_richard_love(q, r, coef, nqubits);
                        cudaq::spin_op number = seeley_richard_love(p, p, 1.0, nqubits);
                        bk_hamiltonian += number * excitation;
                    }
                }
            }
        }
    }

    for (std::size_t p = 0; p < nqubits; p++) {
        for (std::size_t q = 0; q < p; q++) {
            for (std::size_t r = 0; r < q; r++) {
                for (std::size_t s = 0; s < r; s++) {
                    auto coef_pqrs = -two_body_coef(hpqrs, p, q, r, s);
                    if (std::abs(coef_pqrs) > 0.0) {
                        bk_hamiltonian +=  hermitian_one_body_product(p, q, r, s, coef_pqrs, nqubits);
                    }
                    auto coef_prqs = -two_body_coef(hpqrs, p, r, q, s);
                    if (std::abs(coef_prqs) > 0.0) {
                        bk_hamiltonian +=  hermitian_one_body_product(p, r, q, s, coef_prqs, nqubits);
                    }
                    auto coef_psqr = -two_body_coef(hpqrs, p, s, q, r);
                    if (std::abs(coef_psqr) > 0.0) {
                        bk_hamiltonian +=  hermitian_one_body_product(p, s, q, r, coef_psqr, nqubits);
                    }
                }
            }
        }
    }

    bk_hamiltonian += constant_term * cudaq::spin::i(0);
    return bk_hamiltonian;
}
} // namespace cudaqx::solvers
