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
