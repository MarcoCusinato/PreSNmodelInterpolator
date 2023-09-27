import numpy as np

iron_peak_list = ['fe', 'co', 'ni', 'cr', 'mn']
iron_peak_Z = {'fe': 26, 'co': 27, 'ni': 28, 'cr': 24, 'mn': 25}


def find_header_line(lines):
    """
    returns the header line of the species file
    """
    for line in lines:
        if 'nt1' in line.casefold():
            li = line.casefold().split()
            index = li.index('nt1')
            return li[index:]
    raise ValueError('No header line found in species file')

def sum_species(species, file_header, Amin, Amax, remove_element = ''):
    """
    sums up the mass fraction of the species with mass number between Amin and Amax
    """
    Amin = int(abs(Amin))
    Amax = int(abs(Amax))
    assert Amin <= Amax, 'Amin must be smaller than Amax'
    assert len(str(Amin)) == len(str(Amax)), 'Amin and Amax must have the same number of digits'
    out_specie = 0
    if len(str(Amin)) == 1:
        for (index, element) in zip(range(len(file_header)), file_header):
            try:
                int(element[-3:])
            except:
                try:
                    int(element[-2:])
                except:
                    if int(element[-1]) <= Amax and int(element[-1]) >= Amin:
                        out_specie += species[:, index]
    else:
        for (index, element) in zip(range(len(file_header)), file_header):
            try:
                int(element[-3:])
            except:
                try:
                    if int(element[-2:]) <= Amax and int(element[-2:]) >= Amin:
                        if element == remove_element:
                            continue
                        out_specie += species[:, index]
                except:
                    continue
    return out_specie

def sum_iron_species(species, file_header, species_type):
    """
    sums up the mass fraction of the species near the iron peak
    fe54: 2*Z+2 plus fe56
    ni56: < 2*Z + 2
    Fe: > 3+2*Z
    """
    assert species_type in ['fe54', 'ni56', 'Fe'], 'type must be either fe54, ni56 or Fe'
    fe_species = 0
    for (index, element) in zip(range(len(file_header)), file_header):
        if element[:2] not in iron_peak_list:
            continue
        if species_type == 'fe54':
            if 2 * iron_peak_Z[element[:2]] + 2 == int(element[-2:]) or \
                2 * iron_peak_Z[element[:2]] + 3 == int(element[-2:]) or \
                element == 'fe56':
                fe_species += species[:, index]
               
        elif species_type == 'ni56':
            if int(element[-2:]) < 2 + 2 * iron_peak_Z[element[:2]]:
                fe_species += species[:, index]
        elif species_type == 'Fe':
            if int(element[-2:]) > 3 + 2 * iron_peak_Z[element[:2]] and element != 'fe56':
                fe_species += species[:, index]
                
    return fe_species


def convert_species(out_array, species, file_lines):
    """
    Sums up the chemical abundances of the species file and returns them in
    a way consistent with the KEPLER output with reduced species
    """

    header = find_header_line(file_lines)
    sum_species(species, header, 48, 51)
    out_array[:, 0] = species[:, header.index('nt1')]
    out_array[:, 1] = species[:, header.index('h2')]
    out_array[:, 2] = species[:, header.index('he3')]
    out_array[:, 3] = sum_species(species, header, 2, 5)
    out_array[:, 4] = species[:, header.index('c12')]
    out_array[:, 5] = species[:, header.index('n14')]
    out_array[:, 6] = species[:, header.index('o16')]
    out_array[:, 7] = species[:, header.index('ne20')]
    out_array[:, 8] = sum_species(species, header, 23, 28, 'si28')
    out_array[:, 9] = species[:, header.index('si28')]
    out_array[:, 10] = sum_species(species, header, 29, 35)
    out_array[:, 11] = sum_species(species, header, 36, 39)
    out_array[:, 12] = sum_species(species, header, 40, 43)
    out_array[:, 13] = sum_species(species, header, 44, 47)
    out_array[:, 14] = sum_species(species, header, 48, 51)
    out_array[:, 15] = species[:, header.index('fe52')]
    out_array[:, 16] = sum_iron_species(species, header, 'fe54')
    out_array[:, 17] = sum_iron_species(species, header, 'ni56')
    out_array[:, 18] = species[:, header.index('fe56')]
    out_array[:, 19] = sum_iron_species(species, header, 'Fe')
    return out_array