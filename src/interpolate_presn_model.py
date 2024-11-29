import numpy as np
from src.all_species_src.species_conv import convert_species
import json
import os
import platform
import pandas as pd

class unit_converter:
    """
    Class in which the basic conversions are stored.
    """
    def __init__(self):
        self.msol = 1.98847e33
        self.rsol = 6.957e10
    
    def to_m_sol(self, quantity):
        return quantity / self.msol
    
    def to_r_sol(self, quantity):
        return quantity / self.rsol
    
    def to_cm(self, quantity):
        return quantity * self.rsol
    
    def to_g(self, quantity):
        return quantity * self.msol

class InterpolatePresnModel:
    """
    Class InterpolatePresnModel, the purpose of this class is to reshape a pre-supernova model
    in order to be used by Aenus. It also interpolates the quantities in order to have a more
    uniform grid.
    """
    def __init__(self, file_path, models_properties_path, save_path, **kwargs):
        """
        Class initialization, it takes as input the path of the model and the name of the model.
        file_path: str, path to the model file.
        models_properties_path: str, path to the json file containing the properties of the models.
        save_path: str, path where the model will be saved.
        It also takes as input the following keyword arguments:
        - rmin: minimum radius of the grid (default: 0)
        - rmax: maximum radius of the grid (default: 1e13)
        - rmiddle: middle radius of first grid cell (default: 4e4)
        - ngrid: number of grid cells (default: 16000)
        - ftype: KEPLER or MESA (default: None, the format is automatically detected)
        """
        self.u = unit_converter()
        self.rmin = None
        self.rmax = None
        self.rmiddle = None
        self.ngrid = None
        self.ftype = None
        ## SET KWARGS VALUES
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.file_path = file_path
        self.json_file_path = models_properties_path
        ## GET FILE NAME
        self.model_name, self.paper_name = self.__find_file_name()
        ## SET RESULT PATH
        self.result_path = self.__set_result_path(save_path)
        ## LOAD JSON FILE
        self.model_properties = self.__load_properties()
        ## AUTO-DETECT FILE FORMAT
        self.format, self.file_lines, self.header_lines = self.__find_file_format(self.ftype)
        self.footer_lines = self.__find_footer()
        ## DETECT BFIELD PRESENCE
        self.has_bfield = self.__bfield_finder()
        ## LOAD DATA
        self.data = self.__load_data()
        self.thermo, self.nuclei = self.__create_quantities_arrays()
        self.__order_data()
        ## DEFINE THE GRID
        self.__define_grid()
        ## CREATE THE GRID
        self.radius, self.theta = self.__create_grid()
        ## INTERPOLATE THE MODEL
        self.__interpolate_model()
        ## CALCULATE THE MASS
        self.mass = self.__calculate_mass()
        ## SAVE THE MODEL
        self.__save_interpolated_model()
        

    def __set_result_path(self, save_path):
        """
        Method that creates the folder where the model will be saved.
        The folder is located in the results folder and it has the same name of the model.
        """
        print('Creating folder...')
        if save_path is None:
            result_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../results')
            if not os.path.exists(result_path):
                os.mkdir(result_path)
            result_path = os.path.join(result_path, self.model_name)
            if not os.path.exists(result_path):
                os.mkdir(result_path)
        else:
            assert os.path.exists(save_path), 'The path provided does not exist.'
            result_path = os.path.join(save_path, self.model_name)
            if not os.path.exists(result_path):
                os.mkdir(result_path)
        print('Folder created')
        return result_path
 
    def __find_file_name(self):
        """
        Method that finds the name of the model from the file name.
        """
        ## Splitting the path to the file according to the OS##
        if platform.system() == 'Windows':
            path = os.path.abspath(self.file_path).split('\\')
        else:
            path = os.path.abspath(self.file_path).split('/')
        ## Removing the extension and the other not useful stuff##
        extensions = ['mso_final_profile.data', 'mso_final_profile.txt', '.txt',
                      '@presn', '.presn', '.presn_structure', '_presn']
        ## find paper name from folder name ##
        try:
            int(path[-2][-2:])
            paper_name = path[-2]
        except:
            paper_name = input('Please enter the paper name: ')
        file_name = path[-1]
        for extension in extensions:
            file_name = file_name.replace(extension, '')
        ## Check if the file name ia a number ##
        ## If it is, add the first letter of the paper folder at the beginning ##
        try:
            int(file_name)
            if len(file_name) == 1:
                file_name = '0' + file_name
            file_name = paper_name[0] + file_name
        except:
            pass
        file_name = file_name + '_' + paper_name

        return file_name, paper_name

    def __load_properties(self):
        """
        Method that loads the json file containing the properties of the models.
        return: dict, properties of the model.
        """
        print('Loading json file...')
        assert os.path.exists(self.json_file_path), 'The json file does not exist.'
        with open(self.json_file_path) as f:
            json_properties = json.load(f)
        m_name = self.model_name.replace('_' + self.paper_name, '')
    
        if m_name[0] == self.paper_name[0]:
            m_name = m_name[1:]
    
        for model in json_properties:
            if model['name'] == m_name and \
                model['group'] == self.paper_name:
                print('Json file loaded.')
                return model
        raise ValueError('The model is not present in the json file.')
    
    def __find_file_format(self, format=None):
        """
        Find the format of the file. It could be made either with KEPLER or MESA.
        If the format is not specified, it will be automatically detected.
        param format: str, optional (default = None) Can be either 'KEPLER' or 'MESA'.
        return: str, either 'KEPLER_NEW', 'KEPLER_std', 'KEPLER_std_full' or 'MESA',
                list of str, list of lines in the file,
                int, index of the line where the data starts.
        """
        def find_kepler_version(lines):
            ## KEPLER version finder ##
            """
            To distinguish between the three KEPLER formats, we use the following file characteristics:
            1. New KEPLER files have more than 2 commented lines before the data starts.
            2. Old KEPLER files with full nuclear species have 'ni71' in the line before the data starts.
            """
            for lindex in range(1,len(lines)):
                try:
                    float(lines[lindex].split()[1])
                    line_index = lindex
                    break
                except:
                    pass
            if line_index > 2:
                print('New KEPLER format deteted.')
                format = 'KEPLER_NEW'
            else:
                if 'ni71' in lines[line_index - 1]:
                    print('Standard KEPLER format with full atomic species detected.')
                    format = 'KEPLER_std_full'
                else:
                    print('Standard KEPLER format detected.')
                    format = 'KEPLER_std'
            return format, line_index
        
        def check_mesa_file(lines):
            ## MESA file checker ##
            """
            To check if the file is a MESA file, we check if the line with the initial mass is present.
            """
            correct_MESA = False
            for line in lines:
                if 'initial_mass' in line:
                    correct_MESA = True
                    break
            for lindex in range(1,len(lines)):
               if 'logT' in lines[lindex]:
                    line_index = lindex + 1
                    break
            if not correct_MESA:
                raise ValueError('The file is not a MESA file.')
            return line_index
        
        ### Main part of the method ###
        line_index = 0
        if format in ['MESA', 'KEPLER']:
            if format == 'KEPLER':
                with open(self.file_path, 'r') as f:
                    lines = f.readlines()
                format, line_index = find_kepler_version(lines)
            else:
                line_index = check_mesa_file(lines)
        else:
            print('Auto-detecting file format.\nPlease wait...')
            with open(self.file_path, 'r') as f:
                lines = f.readlines()
            if 'VERSION' in lines[0]:
                format, line_index = find_kepler_version(lines)
            else:
                line_index = check_mesa_file(lines)
                print('MESA format detected.')
                format = 'MESA'
        return format, lines, line_index
    
    def __find_footer(self):
        """
        Find the index of the line where the data ends.
        return: int, number on lines to skip from the end of the file.
        """
        findex = 0
        lines = list(reversed(self.file_lines))
        for lindex in range(len(lines)):
            try:
                float(lines[lindex].split()[0].replace(':',''))
                findex = lindex
                break
            except:
                continue
        return findex

    def __bfield_finder(self):
        """
        Checks the header to find weather the file contains magnetic field data.
        return: bool, True if the file contains magnetic field data, False otherwise.
        """
        print('Auto-detecting magnetic field data.\nPlease wait...')
        for line_index in range(self.header_lines):
            if 'b_r' in self.file_lines[line_index].casefold() or \
                'dynamo_log_b_r' in self.file_lines[line_index].casefold():
                print('Magnetic field data detected.')
                return True
        print('No magnetic field data detected.')
        return False

    def __load_data(self):
        """
        Load the data from the file.
        return: numpy array, data from the file.
        """
        print('Loading data from file.\nPlease wait...')
        if self.header_lines > 0:
            print('Skipping {} header lines.'.format(self.header_lines))
        if self.footer_lines == 1:
            print('Last line contains winds, skipping it.')
        elif self.footer_lines > 1:
            print('Skipping {} footer lines.'.format(self.footer_lines))
        data = np.genfromtxt(self.file_path, skip_header=self.header_lines, skip_footer=self.footer_lines)
        print('Data loaded.')
        return data

    def __order_data(self):
        """
        Order the data in the file.
        param data: numpy array, data from the file.
        return: numpy array, ordered data.
        """
        print('Ordering data.\nPlease wait...')
        ## DELETE CELL NUMBER AND STABILITY
        if self.format == 'KEPLER_std_full':
            self.data = np.delete(self.data, (0, 12), 1)
            self.data = np.nan_to_num(self.data)
            print('\treducing full species to 20.\n\tPlease wait...')
            self.__KEPLER_w_species_order_nuclei_thermo()
        elif self.format == 'KEPLER_std':
            self.data = np.delete(self.data, (0, 12, 13), 1)
            self.data = np.nan_to_num(self.data)
            self.__KEPLER_wo_species_order_nuclei_thermo()
        elif self.format == 'KEPLER_NEW':
            self.data = np.delete(self.data, (0, 1), 1)
            self.data = np.nan_to_num(self.data)
            self.__KEPLER_wo_species_order_nuclei_thermo()
        elif self.format == 'MESA':
            self.data = np.delete(self.data, (0), 1)
            self.data = np.flip(self.data, axis=0)
            self.data = np.nan_to_num(self.data)
            self.data = pd.DataFrame(data=self.data,
                                     columns=np.array(self.file_lines[self.header_lines-1].split()[1:],
                                                      dtype=str))
            print('\tInternal energy data not found, calculating it using a perfect gas equation of state.')
            print('\tHeavy iron peak nuclei (\'Fe\' species) not found, setting them to 0.')
            print('\tPlease wait...')
            self.__MESA_order_nuclei_thermo()
        print('Data ordered.')

    def __create_quantities_arrays(self):
        """
        Create the arrays containing the thermodynamic quantities and 
        chemical abundances.
        param data: numpy array, data from the file.
        return: numpy array, thermodynamic quantities,
                numpy array, nuclei.
        """
        nuclei = np.zeros((self.data.shape[0], 21))
        if self.has_bfield:
            thermo = np.zeros((self.data.shape[0], 12))
        else:
            thermo = np.zeros((self.data.shape[0], 10))
        return thermo, nuclei
    
    def __KEPLER_w_species_order_nuclei_thermo(self):
        """
        
        """
        ## Thermodynamic quantities ##
        self.thermo[:, 0] = self.data[:, 1] # radius
        self.thermo[:, 1:3] = self.data[:, 3:5] # density and temperature
        self.thermo[:, 3] = self.data[:, 10] # Ye
        self.thermo[:, 4] = self.data[:, 5] # pressure
        self.thermo[:, 5] = self.data[:, 7] # entropy
        self.thermo[:, 6] = self.data[:, 6] # internal energy
        self.thermo[:, 7] = self.data[:, 9] # Abar
        self.thermo[:, 8] = self.data[:, 2] # velocity
        self.thermo[:, 9] = self.data[:, 8] # omega
        if self.has_bfield:
            self.thermo[:, 10] = self.data[:, -1] # toroidal magnetic field
            self.thermo[:, 11] = self.data[:, -2] # poloidal magnetic field
        ## Nuclei ##
        self.nuclei[:, 0] = self.data[:, 1]
        self.nuclei[:, 1:] = convert_species(self.nuclei[:, 1:], self.data[:, 11:], self.file_lines)

    def __KEPLER_wo_species_order_nuclei_thermo(self):
        """
        Diveìdes the data from the KEPLER standard format in an array containing the
        thermodynamic quantities and one containing the chemical abundances.
        """
    
        ## Thermodynamic quantities ##
        self.thermo[:, 0] = self.data[:, 1] # radius
        self.thermo[:, 1:3] = self.data[:, 3:5] # density and temperature
        self.thermo[:, 3] = self.data[:, 10] # Ye
        self.thermo[:, 4] = self.data[:, 5] # pressure
        self.thermo[:, 5] = self.data[:, 7] # entropy
        self.thermo[:, 6] = self.data[:, 6] # internal energy
        self.thermo[:, 7] = self.data[:, 9] # Abar
        self.thermo[:, 8] = self.data[:, 2] # velocity
        self.thermo[:, 9] = self.data[:, 8] # omega
        if self.has_bfield:
            self.thermo[:, 10] = self.data[:, -1] # toroidal magnetic field
            self.thermo[:, 11] = self.data[:, -2] # poloidal magnetic field
        ## Nuclei ##
        self.nuclei[:, 0] = self.data[:, 1]
        missing_species_indices =self.__check_KEPLER_species()
        if len(missing_species_indices) == 0:
            if self.has_bfield:
                self.nuclei[:, 1:] = self.data[:, 11:-2]
            else:
                self.nuclei[:, 1:] = self.data[:, 11:]
        elif len(missing_species_indices) == 1:
            if self.has_bfield:
                self.nuclei[:, 1:missing_species_indices[0] + 1] = self.data[:, 11:11 + missing_species_indices[0]]
                self.nuclei[:, missing_species_indices[0] + 2:] = self.data[:, 11 + missing_species_indices[0]:-2]
            else:
                self.nuclei[:, 1:missing_species_indices[0] + 1] = self.data[:, 11:11 + missing_species_indices[0]]
                self.nuclei[:, missing_species_indices[0] + 2:] = self.data[:, 11 + missing_species_indices[0]:]
        else:
            raise ValueError('More than one species missing.')
              
    
    def __check_KEPLER_species(self):
        """
        Method that checks the species in the KEPLER file.
        """
        species = ['nt1', 'h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24', 
                   'si28','s32', 'ar36', 'ca40', 'ti44', 'cr48', 'fe52', 'fe54', 'ni56',
                   'fe56', '\'fe\'']
        for line_index in range(self.header_lines):
            if 'h1' in self.file_lines[line_index].casefold():
                if 'neutrons' in self.file_lines[line_index].casefold():
                    species_in_file = self.file_lines[line_index].replace('neutrons', 'nt1').split()
                else:
                    species_in_file = self.file_lines[line_index].split()
                species_in_file = species_in_file[species_in_file.index('nt1'):]
                break
        ## Check if all species are present ##
        missing_elements = list(set(species) - set(species_in_file))
        if len(missing_elements) == 0:
            self.comment = ''
            return []
        else:
            self.comment = 'Missing species: ' + ', '.join(missing_elements)
            missing_species_indices = []
            for el in missing_elements:
                missing_species_indices.append(species.index(el))
            return missing_species_indices
            

    def __MESA_order_nuclei_thermo(self):
        """
        Diveìdes the data from the MESA format in an array containing the
        thermodynamic quantities and one containing the chemical abundances.
        Since no 'Fe' group is provided it is left to 0.
        """
        ## Thermodynamic quantities ##
        self.thermo[:, 0] = self.u.to_cm(np.array(10 ** self.data['logR'])) # radius
        self.thermo[:, 1] = np.array(10 ** self.data['logRho']) # rho
        self.thermo[:, 2] = np.array(10 ** self.data['logT']) # temperature
        self.thermo[:, 3] = np.array(self.data['ye']) # Ye
        self.thermo[:, 4] = np.array(10 ** self.data['logP']) # pressure
        self.thermo[:, 5] = np.array(self.data['entropy']) # entropy
        try:
            self.thermo[:, 6] = np.array(self.data['energy']) # internal energy
        except:
            print('Exception')
            self.thermo[:, 6] = self.thermo[:, 4] / ( np.array(self.data['csound']) ** 2 * \
                                                 self.thermo[:, 1] / self.thermo[:, 4] -1) # internal energy
        self.thermo[:, 7] = np.array(self.data['abar']) # Abar
        self.thermo[:, 8] = np.array(self.data['velocity']) # velocity
        self.thermo[:, 9] = np.array(self.data['omega']) # omega
        if self.has_bfield:
            self.thermo[:, 10] = np.array(10 ** self.data['dynamo_log_B_phi']) # toroidal magnetic field
            self.thermo[:, 11] = np.array(10 ** self.data['dynamo_log_B_r']) # poloidal magnetic field
        ## Nuclei ##
        keys = ['neut', 'h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24',
                'si28', 's32', 'ar36', 'ca40', 'ti44', 'cr48', 'fe52', 'fe54', 'ni56',
                'fe56']
        self.nuclei[:, 0] = self.thermo[:, 0]
        for i, key in enumerate(keys):
            self.nuclei[:, i + 1] = np.array(self.data[key])
        
        
    def __define_grid(self):
        """
        Method that defines the grid of the model, if no grid is provided it uses the default values.
        It also checks if the grid is valid.
        """
        if any([self.rmin is None, self.rmax is None, self.rmiddle is None, self.ngrid is None]):
            print('No grid provided, using default values.')
            self.rmin = 0
            self.rmiddle = 4e4
            if self.nuclei[-1, 0] <= 1e13:
                self.rmax = self.nuclei[-1, 0]
            else:
                self.rmax = 1e13
            self.ngrid = 16000
        assert self.rmin < self.rmiddle < self.rmax, 'rmin < rmiddle < rmax'
        assert self.ngrid > 0, 'ngrid > 0'
        assert self.rmin >= 0, 'rmin >= 0'
    
    def __create_grid(self):
        """
        Method that creates the grid of the model.
        So far only 1D models are supported.
        """
        print('Creating new grid...')
        r = np.zeros((self.ngrid, 4))
        r[:, 0] = np.arange(1, self.ngrid + 1)
        r[:, 3] = np.logspace(np.log10(2 * self.rmiddle - self.rmin), np.log10(self.rmax), self.ngrid)
        r[0, 1] = self.rmin
        r[1:, 1] = r[:-1, 3]
        r[:, 2] = 0.5 * (r[:, 3] + r[:, 1])
        theta = np.array([0, 0.25])
        print('Grid created.')
        return r, theta[None, :]
    
    def __interpolate_model(self):
        """
        Method that interpolates the quantities on the new grid.
        If the old grid does not start from 0, the first two points are interpolated with a parabola.
        """
        print('Interpolating model...')
        self.nuclei_interp = np.zeros((self.ngrid, self.nuclei.shape[1] + 1))
        self.thermo_interp = np.zeros((self.ngrid, self.thermo.shape[1] + 1))
        self.nuclei_interp[:, 0] = np.arange(1, self.ngrid + 1)
        self.nuclei_interp[:, 1] = self.radius[:, 2]
        self.thermo_interp[:, 0] = self.nuclei_interp[:, 0]
        self.thermo_interp[:, 1] = self.radius[:, 2]
        if self.rmin < self.nuclei[0, 0]:
            extrapolation_index = np.argmax(self.radius[:,2] > self.nuclei[0, 0]) + 1
            for i in range(1, self.nuclei.shape[1]):
                a, b, c = self.__parabolic_coefficient(self.nuclei[:, 0], self.nuclei[:, i])
                self.nuclei_interp[:extrapolation_index, i+1] = a * self.radius[:extrapolation_index, 2] ** 2 + b * self.radius[:extrapolation_index, 2] + c
                self.nuclei_interp[extrapolation_index:, i+1] = np.interp(self.radius[extrapolation_index:, 2], self.nuclei[:, 0], self.nuclei[:, i])
            for i in range(1, self.thermo.shape[1]):
                
                a, b, c = self.__parabolic_coefficient(self.thermo[:, 0], self.thermo[:, i])
                self.thermo_interp[:extrapolation_index, i+1] = a * self.radius[:extrapolation_index, 2] ** 2 + b * self.radius[:extrapolation_index, 2] + c
                self.thermo_interp[extrapolation_index:, i+1] = np.interp(self.radius[extrapolation_index:, 2], self.thermo[:, 0], self.thermo[:, i])
        else:
            for i in range(1, self.nuclei.shape[1]):
                self.nuclei_interp[:, i+1] = np.interp(self.radius[:, 2], self.nuclei[:, 0], self.nuclei[:, i])
            for i in range(self.thermo.shape[1]):
                self.thermo_interp[:, i+1] = np.interp(self.radius[:, 2], self.thermo[:, 0], self.thermo[:, i])
        if self.has_bfield:
            self.thermo_interp[:, -1] = np.where(self.thermo_interp[:, -1] < 0, 0.0, self.thermo_interp[:, -1])
            self.thermo_interp[:, -2] = np.where(self.thermo_interp[:, -2] < 0, 0.0, self.thermo_interp[:, -2])
        print('Model interpolated')
    
    def __calculate_mass(self):
        """
        Method that calculates the mass of the model, encloded in the grid.
        """
        print('Calculating mass...')
        dr = 4 * np.pi * (self.radius[:, 3] ** 3 - self.radius[:, 1] ** 3 ) / 3
        mass = np.sum(dr * self.thermo_interp[:, 2]) / 1.988e33
        print('Mass calculated.')
        return mass
    
    def __create_nuclei_text(self):
        """
        Method that creates the text to be written in the nuclei.pars file.
        """
        if self.format == 'MESA':
            with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nuclei_file/nuclei_MESA.pars')) as f:
                nuclei_text = f.readlines()
            nuclei_text[15] = nuclei_text[15].replace('MESA', self.paper_name)
        else:
            with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nuclei_file/nuclei_Heger.pars')) as f:
                nuclei_text = f.readlines()
            nuclei_text[15] = nuclei_text[15].replace('WHW2002', self.paper_name)
        return nuclei_text

    def __update_model_list(self):
        """
        Method that updates the model list in the json file.
        """
        ## custom sort function ##
        if "ZAMS_mass1" in self.model_properties.keys():
            sorted_keys = ["name", "group", "original_model", "ZAMS_mass1",
                    "ZAMS_mass2", "enclosed_mass", "mass", "xi15", "xi175",
                    "xi25", "star_type", "metallicity", "omg", "btor", "bpol",
                    "comment"]
        else:
            sorted_keys = ["name", "group", "original_model", "ZAMS_mass",
                    "enclosed_mass", "mass", "xi15", "xi175", "xi25", 
                    "star_type", "metallicity", "omg", "btor", "bpol",
                    "comment"]
        ## Add key to the json file ##
        self.model_properties['original_model'] = self.model_properties['name']
        self.model_properties['name'] = self.model_name
        self.model_properties['enclosed_mass'] = '{:.1f}'.format(self.mass)
        if self.format == 'MESA':
            self.model_properties['comment'] += ' Internal energy calculated from a perfect gas equation of state. Only 19 element species are provided.'
        if self.comment != '':
            self.model_properties['comment'] += ' ' + self.comment + '.'
        ## sort the dictionary ##
        self.model_properties = {k: self.model_properties[k] for k in sorted_keys}
        json_path = self.result_path.replace(self.model_name, 'models_list.json')
        if os.path.exists(json_path):
            with open(json_path, 'r') as f:
                json_data = json.load(f)
            json_data.append(self.model_properties)
        else:
            json_data = [self.model_properties]
            print('Json file not found, creating a new one...')
        print('Updating json file...')

        ## Clear all duplicates in the json ##
        result = list()
        items_set = set()

        for js in json_data:
            # only add unseen items (referring to 'title' as key)
            if not js['name'] in items_set:
                # mark as seen
                items_set.add(js['name'])
                # add to results
                result.append(js)
        json_data = result
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=4)
        
    
    def __save_interpolated_model(self):
        """
        Method that saves the interpolated model.
        It saves the nuclei, the thermodynamic quantities and the grid in four different files.
        It also produces two additional files, one containing the thermodynamic quantities information
        and one containing the chemical abundances information.
        """
        print('Saving model...')
        print('Saving nuclei...')
        nuclei_fmt = '\t%d'
        for i in range(self.nuclei_interp.shape[1] - 1):
            nuclei_fmt += '\t%.20E'
        np.savetxt(os.path.join(self.result_path, 'nuclei.dat'), self.nuclei_interp, fmt=nuclei_fmt)
        with open(os.path.join(self.result_path, 'nuclei.pars'), 'w') as f:
            f.writelines(self.__create_nuclei_text())
        
        print('Saving thermodynamic quantities...')
        thermo_fmt = '\t%d'
        for i in range(self.thermo_interp.shape[1] - 1):
            thermo_fmt += '\t%.20E'
        np.savetxt(os.path.join(self.result_path, 'star.dat'), self.thermo_interp, fmt=thermo_fmt)
        with open(os.path.join(self.result_path, 'star.txt'), 'w') as f:
            f.write(self.__produce_text())
        with open(os.path.join(self.result_path, 'Heger.pars'), 'w') as f:
            f.write(self.__produce_Heger_pars())
        print('Saving grid...')
        np.savetxt(os.path.join(self.result_path, 'initial_model.x.dat'), self.radius, fmt = '\t%d\t%.20E\t%.20E\t%.20E')
        np.savetxt(os.path.join(self.result_path, 'initial_model.y.dat'), self.theta, fmt = '\t%d\t%.20E')
        self.__update_model_list()
        print('Model saved.')

    def __produce_text(self):
        """
        Method that produces the text to be written in the star.txt file.
        It saves the minimum and maximum of the gird, the number of grid cells and the mass of the model.
        It also provides the colums of the data file.
        """
        output_text = '----------------------------------------------------------\n' + \
                    ' Heger model no. ' + self.model_name + '\n' + \
                    '----------------------------------------------------------\n' + \
                    ' grid of ngrid zones ranging from rmin to rmax\n' + \
                    ' rmin = ' + str(self.rmin) + '\n' + \
                    ' rmax = ' + str(self.rmax) + '\n' + \
                    ' ngrid = ' + str(self.ngrid) + '\n' + \
                    ' mass = ' + str(self.mass) + '\n' + \
                    '----------------------------------------------------------\n' + \
                    ' data written to file ' + self.model_name +  '.dat are ::\n'
        if self.has_bfield:
            output_text += '   rho, tem, y_e, pre, ent, erg, abr, vel, omg, btor, bpol\n' + \
                        '----------------------------------------------------------\n'
        else:
            output_text += '   rho, tem, y_e, pre, ent, erg, abr, vel, omg\n' + \
                        '----------------------------------------------------------\n'
        if self.format == 'MESA':
            output_text += '   Warning: internal energy derived with perfect gas\n' + \
                        'EOS: U = P / (gamma - 1)' + \
                        '----------------------------------------------------------\n'
        return output_text
    
    def __produce_Heger_pars(self):
        """
        Method that produces the namelist to be written in the Heger.pars file.
        """
        str_data = '\t1\t' + str(self.ngrid) + '\n' + \
                   '\t1\t1\n' + \
                   '\t1\t1\n' + \
                   '&Heger_pars\n' + \
                   '  heger_nx =\t' + str(self.ngrid) + '\n' + \
                   '/'
        return str_data 

    ## UTILITY METHODS
    def __parabolic_coefficient(self, x, y):
        """
        Method that calculates the coefficients of a parabola given three points.
        """
        x_points = x[:4]
        y_points = y[:4]
        return np.polyfit(x_points, y_points, 2)
