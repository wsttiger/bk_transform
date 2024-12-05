#include <set>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iterator>

#include "spin_op.h"
#include "tensor.h"

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

    cudaq::spin_op seeley_richard_love_result = 0;
    // Case 0
    if (i == j) {
        std::cout << "Case 0\n";
        cudaq::spin_op ops;
        for (int index : occupation_set(i)) {
            ops *= cudaq::spin::z(index);
        }
        seeley_richard_love_result +=  coef * 2.0 * cudaq::spin_op();
        seeley_richard_love_result += -coef * 2.0 * ops;
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
            c0 = -imag_i * coef;
        } else {
            c0 =          -coef;
            c1 = -imag_i * coef;
            c1 =  imag_i * coef;
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

cudaq::spin_op generate(const double constant,
                        const cudaqx::tensor<> &hpq,
                        const cudaqx::tensor<> &hpqrs) {

    auto nqubits = hpq.shape()[0];
    cudaq::spin_op bk_hamiltonian = 0;
    for (std::size_t p = 0; p < nqubits; p++) {
        if (std::abs(hpq.at({p,p})) > 0.0) {
            bk_hamiltonian += seeley_richard_love(p, p, hpq.at({p,p}), nqubits); 
        }
        for (std::size_t q = 0; q < p; q++) {
            if (std::abs(hpq.at({p,q})) > 0.0) {
                bk_hamiltonian += seeley_richard_love(p, q, std::conj(hpq.at({p,q})), nqubits); 
            }
        }
    }
    return bk_hamiltonian;
}

int main() {
    int n_qubits = 6;
    cudaq::spin_op op;
    cudaq::spin_op identity;
    op += cudaq::spin::x(4)*cudaq::spin::z(5);
    op -= identity;
    {
        n_qubits = 20;
        for (std::size_t i = 0; i < n_qubits; i++) {
            for (std::size_t j = 0; j < n_qubits; j++) {
                auto result = seeley_richard_love(i, j, 4.0, n_qubits);
                std::cout << "i: " << i << "  j: " << j << "\n";
                result.dump();
     
            }
        }
    }
#if 0
    {
        auto result = seeley_richard_love(2, 2, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(1, 2, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(1, 3, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(1, 4, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(3, 4, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(3, 5, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(4, 5, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(2, 3, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(2, 4, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(2, 5, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
    {
        auto result = seeley_richard_love(2, 1, 4.0, n_qubits);
        result.dump();
        std::cout << "\n";
    }
#endif
    return 0;
}

