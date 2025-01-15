# File: urqrd_functions.pyx

import numpy as np
cimport numpy as np
from libc.math cimport sqrt
from libc.stdlib cimport malloc, free
from cpython.list cimport  PyList_Append, PyList_GET_SIZE

def urqrd_main(np.ndarray[np.complex128_t, ndim=1] data, int rank, int iterations):
    cdef int order = int(data.size/2)
    cdef int N = len(data) - order + 1
    cdef int i
    cdef np.ndarray[np.float64_t, ndim=2] random_numbers
    cdef np.ndarray[np.complex128_t, ndim=2] data_multiplied_by_random_q
    cdef np.ndarray[np.complex128_t, ndim=2] first_order_q
    
    for i in range(iterations):
        random_numbers = np.random.normal(size=(N, rank))
        data_multiplied_by_random_q = urqrd_core(data, random_numbers)
        first_order_q, _ = np.linalg.qr(FastHankel_prod_mat_mat(data, random_numbers))
        data = FastHankel_2dt(data, first_order_q, data_multiplied_by_random_q)
    
    if data.dtype == "float":
        data = np.real(data)
    
    return data

cdef np.ndarray[np.complex128_t, ndim=2] urqrd_core(np.ndarray[np.complex128_t, ndim=1] data, np.ndarray[np.float64_t, ndim=2] random_numbers):
    cdef np.ndarray[np.complex128_t, ndim=2] data_multiplied_by_random_numbers = FastHankel_prod_mat_mat(data, random_numbers)
    cdef np.ndarray[np.complex128_t, ndim=2] first_order_q
    cdef np.ndarray[np.complex128_t, ndim=2] r
    first_order_q, r = np.linalg.qr(data_multiplied_by_random_numbers)
    cdef np.ndarray[np.complex128_t, ndim=2] data_multiplied_by_random_q = FastHankel_prod_mat_mat(data.conj(), first_order_q).conj().T
    return data_multiplied_by_random_q

cdef np.ndarray[np.complex128_t, ndim=2] FastHankel_prod_mat_mat(np.ndarray[np.complex128_t, ndim=1] data, np.ndarray[np.float64_t, ndim=2] random_numbers):
    cdef int N = random_numbers.shape[0]
    cdef int rank = random_numbers.shape[1]
    cdef int length_of_data = len(data)
    cdef int N_star = length_of_data - N + 1
    cdef np.ndarray[np.complex128_t, ndim=2] d = np.zeros((N_star, rank), dtype=complex)
    cdef int r
    cdef np.ndarray[np.float64_t, ndim=1] column_of_random_nums
    
    for r in range(rank):
        column_of_random_nums = random_numbers[:, r]
        d[:, r] = FastHankel_prod_mat_vec(column_of_random_nums, data)
    return d

cdef np.ndarray[np.complex128_t, ndim=1] FastHankel_prod_mat_vec(np.ndarray[np.float64_t, ndim=1] column_of_random_numbers, np.ndarray[np.complex128_t, ndim=1] data):
    cdef np.ndarray[np.float64_t, ndim=1] column_of_random_numbers_zerofilled = np.concatenate((np.zeros(len(data) - len(column_of_random_numbers)), column_of_random_numbers[::-1]))
    cdef np.ndarray[np.complex128_t, ndim=1] data_fft = np.fft.fft(data)
    cdef np.ndarray[np.complex128_t, ndim=1] random_fft = np.fft.fft(column_of_random_numbers_zerofilled)
    cdef np.ndarray[np.complex128_t, ndim=1] product_of_two_ffts = data_fft * random_fft
    cdef np.ndarray[np.complex128_t, ndim=1] resultant_signal = np.fft.ifft(product_of_two_ffts)
    return np.roll(resultant_signal, +1)[:(len(data) - len(column_of_random_numbers) + 1)]

cdef np.ndarray[np.complex128_t, ndim=1] FastHankel_2dt(np.ndarray[np.complex128_t, ndim=1] data, np.ndarray[np.complex128_t, ndim=2] first_order_q, np.ndarray[np.complex128_t, ndim=2] data_multiplied_by_random_q):
    cdef int data_len_half = first_order_q.shape[0]
    cdef int rank = first_order_q.shape[1]
    cdef int data_len_half_plus1 = data_multiplied_by_random_q.shape[1]
    cdef np.ndarray[np.complex128_t, ndim=1] sum_of_all_data = np.zeros(len(data), dtype=complex)
    cdef int r
    cdef np.ndarray[np.complex128_t, ndim=1] row
    cdef np.ndarray[np.complex128_t, ndim=1] zerofilled_first_order_q
    cdef np.ndarray[np.complex128_t, ndim=1] zerofilled_first_order_q_multi_data
    
    for r in range(rank):
        row = data_multiplied_by_random_q[r, :]
        zerofilled_first_order_q = np.concatenate((np.zeros(data_len_half_plus1-1), first_order_q[:, r], np.zeros(data_len_half_plus1-1)))
        zerofilled_first_order_q_multi_data = FastHankel_prod_mat_vec(row[::-1], zerofilled_first_order_q)
        sum_of_all_data += zerofilled_first_order_q_multi_data
    return sum_of_all_data * vec_mean(data_len_half, len(data))

cdef np.ndarray[np.float64_t, ndim=1] vec_mean(int M, int L):
    cdef np.ndarray[np.float64_t, ndim=1] vec_prod_diag = np.array([1.0 / (i + 1) for i in range(M)], dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] vec_prod_middle = np.full(L - 2*M, 1.0 / M, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] vec_mean_prod_tot = np.concatenate((vec_prod_diag, vec_prod_middle, vec_prod_diag[::-1]))
    return vec_mean_prod_tot