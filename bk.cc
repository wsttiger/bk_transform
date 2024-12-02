#include <set>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iterator>

#include "spin_op.h"

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

/**
 * Algebraic expressions for general products of the form a_i^d a_j term in
 * the Bravyi-Kitaev basis. These expressions vary in form depending on the
 * parity of the indices i and j, as well as on the overlaps between the parity
 * and update sets of the indices
 */
void seeley_richard_love(int i, int j, double coef, int n_qubits) {

    using double_complex = std::complex<double>;

    coef *= 0.25;

    cudaq::spin_op seeley_richard_love_result = 0.0 * cudaq::spin::i(0);
    // Case 0
    if (i == j) {  // Simplifies to the number operator
        cudaq::spin_op ops;
        // Create occupation set operator
        for (int index : occupation_set(i)) {
            ops *= cudaq::spin::z(index);
        }
        seeley_richard_love_result +=  coef * 2.0 * cudaq::spin_op();
        seeley_richard_love_result += -coef * 2.0 * ops;
    }

    // Case 1
    else if (i % 2 == 0 and j % 2 == 0) {
        // Create padding operators
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

        // Combine padding operators
        cudaq::spin_op left_pad = x_pad * y_pad * z_pad;

        // Create four operator combinations
        cudaq::spin_op op1 = left_pad * cudaq::spin::y(j) * cudaq::spin::x(i);
        cudaq::spin_op op2 = left_pad * cudaq::spin::x(j) * cudaq::spin::y(i);
        cudaq::spin_op op3 = left_pad * cudaq::spin::x(j) * cudaq::spin::x(i);
        cudaq::spin_op op4 = left_pad * cudaq::spin::y(j) * cudaq::spin::y(i);

        if (i < j) {
            seeley_richard_love_result += double_complex( coef, 0.0) * op1;
            seeley_richard_love_result += double_complex(-coef, 0.0) * op2;
            seeley_richard_love_result += double_complex(0.0, -coef) * op3;
            seeley_richard_love_result += double_complex(0.0, -coef) * op4;
        } else {  // whenever i < j introduce phase of -j
            seeley_richard_love_result += double_complex(0.0, -coef) * op1;
            seeley_richard_love_result += double_complex(0.0,  coef) * op2;
            seeley_richard_love_result += double_complex(-coef, 0.0) * op3;
            seeley_richard_love_result += double_complex(-coef, 0.0) * op4;
        }
    }

    // Case 2
    else if (i % 2 == 1 and j % 2 == 0 and not parity_set(j).contains(i)) {
        cudaq::spin_op x_pad;
        for (int index : U_diff_a_set(i, j, n_qubits)) {
            x_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op y_pad;
        for (int index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }

        // Combine padding operators
        cudaq::spin_op left_pad = x_pad * y_pad;

        cudaq::spin_op right_pad_1;
        std::set<std::size_t> P0_minus_alpha;
        std::set_difference(P0_ij_set(i, j).begin(), 
                            P0_ij_set(i, j).end(), 
                            alpha_set(i, j, n_qubits).begin(), 
                            alpha_set(i, j, n_qubits).end(), 
                            std::inserter(P0_minus_alpha, P0_minus_alpha.begin()));

        for (int index : P0_minus_alpha) {
            right_pad_1 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_2;
        std::set<std::size_t> P2_minus_alpha;
        std::set_difference(P2_ij_set(i, j).begin(), 
                            P2_ij_set(i, j).end(), 
                            alpha_set(i, j, n_qubits).begin(), 
                            alpha_set(i, j, n_qubits).end(), 
                            std::inserter(P2_minus_alpha, P2_minus_alpha.begin()));
        for (int index : P2_minus_alpha) {
            right_pad_2 *= cudaq::spin::z(index);
        }

        std::cout << "[";
        for (auto p : P0_ij_set(i,j)) {
            std::cout << p << " "; 
        }
        std::cout << "]\n\n";
        std::cout << "[";
        for (auto p : P2_ij_set(i,j)) {
            std::cout << p << " "; 
        }
        std::cout << "]\n\n";

        right_pad_1.dump();
        right_pad_2.dump();

        seeley_richard_love_result = left_pad * cudaq::spin::y(j) * cudaq::spin::x(i) * right_pad_1;
        seeley_richard_love_result = left_pad * cudaq::spin::x(j) * cudaq::spin::x(i) * right_pad_1;
        seeley_richard_love_result = left_pad * cudaq::spin::x(j) * cudaq::spin::y(i) * right_pad_2;
        seeley_richard_love_result = left_pad * cudaq::spin::y(j) * cudaq::spin::y(i) * right_pad_2;

if (i < j) {

} else {
}
seeley_richard_love_result.dump();
}
}

int main() {
    int n_qubits = 4;
    cudaq::spin_op op;
    cudaq::spin_op identity;
    op += cudaq::spin::x(4)*cudaq::spin::z(5);
    op -= identity;
    // seeley_richard_love(1, 1, 2.5, n_qubits);
    seeley_richard_love(3, 2, 2.5, n_qubits);

    return 0;
}

