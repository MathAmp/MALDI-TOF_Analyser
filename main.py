import os
import itertools
from typing import List, Tuple, Dict

max_mass = 3500
h2o_mass = 18.015

mass_table: Dict[str, float] = {
    'G': 75.00, 'A': 89.04, 'S': 105.03, 'P': 115.05, 'V': 117.07,
    'T': 119.05, 'C': 121.01, 'L': 131.08, 'I': 131.08, 'N': 132.04,
    'D': 133.03, 'Q': 146.06, 'K': 146.09, 'E': 147.04, 'M': 149.04,
    'H': 155.06, 'F': 165.07, 'R': 174.10, 'Y': 181.06, 'W': 204.08
}

# constant for chem report
pep_dir = 'pep'
alpha = 'alpha.txt'
beta = 'beta.txt'
bsa = 'bsa.txt'
lyso = 'lyso.txt'
oval = 'oval.txt'

# path constant for chem report
alpha_path = os.path.join(pep_dir, alpha)
beta_path = os.path.join(pep_dir, beta)
bsa_path = os.path.join(pep_dir, bsa)
lyso_path = os.path.join(pep_dir, lyso)
oval_path = os.path.join(pep_dir, oval)


def read_file(name: str) -> List[str]:
    """
    read regular file
    :param name: file name
    :return: list of line(str)
    """
    ret = list()
    with open(name, 'r') as file:
        while line := file.readline().strip():
            ret.append(line)
    return ret

# protein peptide sequence
alpha_seq = read_file(alpha_path)[0]
beta_seq = read_file(beta_path)[0]
bsa_seq = read_file(bsa_path)[0]
lyso_seq = read_file(lyso_path)[0]
oval_seq = read_file(oval_path)[0]

# hard coded peak data
a_peak = [1439.791, 1479.667, 1567.782, 1589.722, 2044.853]
b_peak = [1833.711, 1274.734, 2089.932, 2882.583, 1477.706]
c_peak = [874.428, 1045.555, 1428.569, 1675.743, 1804.753]
d_peak = [1858.972, 1687.918, 1581.833, 2283.792, 2788.672]


def get_slice_fragment(sequence: str) -> List[str]:
    """
    ALL possible slice fragment according to trypsin slicing algorithm
    :param sequence: protein peptide sequence
    :return: list of peptide fragment
    """
    slice_position = get_slice_position(sequence)
    return [sequence[first: second] for first, second in itertools.combinations(slice_position, 2)]


def get_ideal_slice_fragment(sequence: str) -> List[str]:
    """
    IDEAL slice fragment according to trypsing slicing algorithm
    :param sequence: protein peptide sequence
    :return: list of peptide fragment
    """
    slice_position = get_slice_position(sequence)
    return [sequence[first: second] for first, second in zip(slice_position, slice_position[1:])]


# just for assignment printing
def assignment_len_limit(seq):
    return 5 <= len(seq) <= 20


# just for assignment printing
def print_fragment(sequence):
    ideal = list(filter(assignment_len_limit, get_ideal_slice_fragment(sequence)))
    total = list(filter(assignment_len_limit, get_slice_fragment(sequence)))
    ideal_tuple = [(seq, get_mass(seq)) for seq in ideal]
    total_tuple = [(seq, get_mass(seq)) for seq in total]
    diff_tuple = set(total_tuple).difference(set(ideal_tuple))
    return list(sorted([(name, round(num, 2)) for name, num in ideal_tuple]))
    #print(list(sorted(diff_tuple)))


def get_slice_position(sequence: str) -> List[int]:
    """
    trypsing algorithm slicing spot. return next value + head + tail
    :param sequence:
    :return:
    """
    ret = set()
    for idx, char in enumerate(sequence):
        if char == 'R' or char == 'K':
            ret.add(idx + 1)
    return sorted({0, len(sequence)} | ret)


def get_mass(sequence: str) -> float:
    """
    return mass of peptide fragment
    :param sequence: peptide fragment code
    :return: mass(unit in Da)
    """
    return sum([mass_table[char] for char in sequence], float()) - h2o_mass * (len(sequence) - 1)


# for console
def make_mass_sequence_map(sequence: str, save_dict: dict) -> None:
    save_dict[get_mass(sequence)] = sequence


# for console
def make_map(all_sequence: List[str], save_dict: dict) -> None:
    for sequence in all_sequence:
        make_mass_sequence_map(sequence, save_dict)


# for convenience of analysis (length limit)
def limit_mass(mass):
    return mass < max_mass


def get_possible_mass(sequence: str) -> List[float]:
    """
    return (total) possible fragment mass list
    :param sequence: protein peptide code
    :return: list of possible mass
    """
    return sorted(filter(limit_mass, map(get_mass, get_slice_fragment(sequence))))


def get_min_abs_distance(mass: float, mass_sequence: List[float]) -> float:
    """
    return min error distance about given mass and mass_sequence(probably theoretical)
    :param mass: experiment single value
    :param mass_sequence: mass sequence (probably theoretical)
    :return:
    """
    return min([abs(mass - x) / x for x in mass_sequence])


def get_min_abs_value(mass: float, mass_sequence: List[float]) -> float:
    """
    return the value that gives min error distance about given mass and mass_sequence(probably theoretical)
    :param mass: experiment single value
    :param mass_sequence: mass sequence (probably theoretical)
    :return:
    """
    return min([(abs(mass - x) / x, x) for x in mass_sequence])[1]


def get_distance(detected_list: List[float], mass_list: List[float]) -> float:
    """
    square sum of min abs distance in detected list
    :param detected_list: experiment value list
    :param mass_list: mass sequence (probably theoretical)
    :return:
    """
    return sum([pow(get_min_abs_distance(d_mass, mass_list), 2) for d_mass in detected_list])


def get_mean_distance(detected_list: List[float], mass_list: List[float]) -> float:
    """
    detected_list X theorical mass list -> square mean of minimum experiment error
    :param detected_list:
    :param mass_list:
    :return:
    """
    return get_distance(detected_list, mass_list) / len(detected_list)


def get_distance_from_sequence(detected_list: List[float], sequence: str) -> float:
    """
    detected_list X sequence(str) -> square mean of minimum experiment error
    :param detected_list: experiment data
    :param sequence:
    :return:
    """
    return get_mean_distance(detected_list, get_possible_mass(sequence))


def shift_detected_mass(detected_list: List[float]) -> List[float]:
    """
    shift float data list - h2o mass
    :param detected_list:
    :return:
    """
    return [x - h2o_mass for x in detected_list]


def read_data(name: str) -> List[float]:
    """
    read data and return float data list
    :param name: file name
    :return: float data list
    """
    return list(map(float, read_file(f"{name}.txt")))


def get_all_exp_data(name: str) -> List[float]:
    """
    read exp data and shift (-1.008)
    :param name: file name
    :return: shifted float data list
    """
    return shift_detected_mass(read_data(name))


alpha_mass = get_possible_mass(alpha_seq)
beta_mass = get_possible_mass(beta_seq)
hem_mass = sorted(set(alpha_mass) | set(beta_mass))
bsa_mass = get_possible_mass(bsa_seq)
lyso_mass = get_possible_mass(lyso_seq)
oval_mass = get_possible_mass(oval_seq)
