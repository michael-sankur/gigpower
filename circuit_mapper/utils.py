from typing import List
import numpy as np


def parse_phases(phase_char_lst: List[str]) -> List[bool]:
    """
    helper function to return a list of phase characters into a boolean triple
    ex: ['1', '3'] -> [True, False True]
    ex: ['a','b'] -> [True, True, False]
    """
    phase_list = [False, False, False]
    for p in phase_char_lst:
        phase_list[get_phase_idx(p)] = True
    return phase_list


def get_phase_idx(phase_char: str) -> int:
    """
    helper function to turn a phase letter into an index, where 'a' = 0
    """
    if phase_char in ['a', 'b', 'c']:
        return ord(phase_char.lower()) - ord('a')
    elif phase_char in ['1', '2', '3']:
        return int(phase_char) - 1
    else:
        raise ValueError(f'Invalid argument for get_phase_idx {phase_char}')


def set_zip_values(dss, zipv: List):
    """sets custom zip values in dss by setting the dss.Loads.zipv() array."""
    # array mapping: [a_z_p, a_i_p, a_pq_p, a_z_q, a_i_q, a_pq_q, min votlage pu]
    for load_name in dss.Loads.AllNames():
        dss.Loads.Name(load_name)
        dss.Loads.Model(8)
        dss.Loads.ZipV(zipv)
        dss.Loads.ZipV()


def pad_phases(matrix: np.ndarray, shape: tuple, phases: List[bool]) -> np.ndarray:
    """
    Helper method to reshape input matrix and set values set to 0
    for phases set to FALSE in phases tuple.
    Input:
        matrix: an nd array
        shape: a 2-element tuple indicating the output matrix's shape
        phases: a tuple of booleans corresponding to phases to set on this matrix (A: T/F, B: T/F, C: T/F)
    Output:
        matrix reshaped with original values, and
        with 0's for all row/column indices corresponding to phases set to FALSE
    """
    # make the return matrix matrix
    ret_mat = np.zeros(shape, dtype=complex)
    vals = iter(matrix.flatten())
    for out_idx in range(shape[0]):
        if len(shape) == 2:
            for col_idx in range(shape[1]):
                if phases[out_idx] and phases[col_idx]:
                    try:
                        ret_mat[out_idx][col_idx] = next(vals)
                    except StopIteration:
                        ("Cannot pad matrix.")
        else:
            try:
                if phases[out_idx]:
                    ret_mat[out_idx] = next(vals)
            except StopIteration:
                ("Cannot pad matrix.")
    return ret_mat


def mask_phases(matrix: np.ndarray, shape: tuple, phases: List[bool]) -> np.ndarray:
    """
    Zeroes out values in input matrix for phases set to FALSE in the phases tuple.
    Input:
        matrix: a 3x3 ndarray
        phases: a tuple of booleans corresponding to phases to set on this matrix (A: T/F, B: T/F, C: T/F)
    Output:
        input matrix with 0's for all row/column indices corresponding to phases set to FALSE
    """
    phase_matrix = np.zeros(shape, dtype=complex)
    for out_idx in range(shape[0]):
        if len(shape) == 2:
            for col_idx in range(shape[1]):
                if phases[out_idx] and phases[col_idx]:
                    phase_matrix[out_idx][col_idx] = 1
        else:
            if phases[out_idx]:
                phase_matrix[out_idx] = 1
    masked = np.multiply(matrix, phase_matrix)
    # change all NaN's to 0
    return np.nan_to_num(masked)
