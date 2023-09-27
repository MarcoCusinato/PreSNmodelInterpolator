from src.interpolate_presn_model import InterpolatePresnModel
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--model-path', type=str, required=True, help='Path to the model file')
parser.add_argument('--json-models-properties-path', type=str, 
                    default='../original/properties.json', 
                    help='Path tho the json file containing the models properties')
parser.add_argument('--save-path', type=str,
                    default='../presn_models', 
                    help='Path to save the interpolated model. If None, the model is saved in the \"reusults\" folder')
parser.add_argument('--rmin', type=float, default=None, help='Minimum of the radius, default is 0.0')
parser.add_argument('--rmax', type=float, default=None, help='Maximum of the radius, default is 1e13')
parser.add_argument('--rmiddle', type=float, default=None, help='Middle radius of first grid cell, default is 4e4')
parser.add_argument('--ngrid', type=int, default=None, help='Number of grid cells, default is 16000')
parser.add_argument('--file-type', type=int, default=None, help='1 (old Heger file format) thermodynamics quantities' + \
                    ' followed by mass fraction of all the elements, 2 (new Heger file format) thermodynamics quantities followed by 20 atomic species')
parser.add_argument('--has-bfield', action='store_true', help='Magnetic field is included in the model')

args = parser.parse_args()
InterpolatePresnModel(file_path=args.model_path,
                      models_properties_path=args.json_models_properties_path,
                      save_path=args.save_path,
                      rmin = args.rmin,
                      rmax = args.rmax,
                      rmiddle = args.rmiddle,
                      ngrid = args.ngrid,
                      ftype = args.file_type)
