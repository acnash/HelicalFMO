# Main entry into HelicalFMO
import argparse
import os
import yaml
from typing import List

from logger_config import get_logger
from controllers.contact_controller import ContactController
from src.controllers.cap_controller import CapController
from src.controllers.fmo_controller import FMOController
from src.controllers.gen_helix_controller import GenHelixController
from src.controllers.rotation_controller import RotationController


def main() -> None:
    logger = get_logger(__name__)
    section_to_run = ["contact_distance", "fmo", "cap", "generate", "rotation"]

    temp_folder = "../temp"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    parser = argparse.ArgumentParser(description="Process TM domain PDB files or folder of files.")

    #parser.add_argument("--file", type=str, required=False,
    #                    help="Path to target file for contact distance and FMO input file generator.")
    #parser.add_argument("--folder", type=str, required=False,
    #                    help="Path to target folder for contact distance and FMO input file generator.")
    #parser.add_argument("--mode", type=str, choices=["contact_distance", "fmo", "cap", "generate", "rotation"],
    #                    required=True, help="Select a mode: contact_distance, fmo, cap, generate, rotation.")
    #parser.add_argument("--output_folder", type=str, required=False)
    #parser.add_argument("--distance_cutoff", type=float, required=False, default=8)
    #parser.add_argument("--ignore_num_start_res", type=int, required=False, default=0)
    #parser.add_argument("--ignore_num_end_res", type=int, required=False, default=0)
    #parser.add_argument("--renum_chains", nargs="+",
    #                    help="List of chain:resid pairs (e.g., A:4 B:99 C:27). Renumbering starts from the resid.")
    #parser.add_argument("--basis", type=str, choices=["STO-3G", "6-31G*"], default="STO-3G", required=False,
    #                    help="FMO basis set: STO-3G (default), 6-31G*")
    #parser.add_argument("--theory", type=str, choices=["HF", "MP2"], default="HF", required=False,
    #                    help="FMO level of theory: HF (default), MP2")
    #parser.add_argument("--seq_a", type=str, required=False, help="First helix sequence (if provided alone, forms a homodimer).")
    #parser.add_argument("--seq_b", type=str, required=False,
    #                    help="Second helix sequence to form a heterodimer")
    #parser.add_argument("--output_file", type=str, required=False, help="File path and file name to save the generated hetero/homodimer as a PDB.")
    #parser.add_argument("--rotation_angle", type=int, default=20, required=False, help="The rotation angle. Default is 20 degrees. ")

    parser.add_argument("--input", type=str, required=True, help="YAML input file.")
    args = parser.parse_args()
    input_file = args.input

    with open(input_file, "r") as file:
        yaml_dict = yaml.safe_load(file)

    for section in section_to_run:
        if section in yaml_dict:
            config_section = yaml_dict[section]

            #if <mode X>
                #generate the controller
                #pass config_section into the new controller
                #run the controller

    #file_location = args.file
    #folder_location = args.folder
    #output_folder = args.output_folder
    #distance_cutoff = args.distance_cutoff
    #ignore_num_start_res = args.ignore_num_start_res
    #ignore_num_end_res = args.ignore_num_end_res
    #mode = args.mode
    #renum_chains_list = args.renum_chains
    #basis = args.basis
    #theory = args.theory
    #sequence_A = args.seq_a
    #sequence_B = args.seq_b
    #output_file = args.output_file
    #rotation_angle = args.rotation_angle

    #if ignore_num_start_res < 0:
    #    print(f"Error: --ignore_num_start_res must be >= 0")
    #    logger.error(f"Error: --ignore_num_start_res must be >= 0")
    #    return

    #if ignore_num_end_res < 0:
    #    print(f"Error: --ignore_num_end_res must be >= 0")
    #    logger.error(f"Error: --ignore_num_end_res must be >= 0")
    #    return

    # these should be pushed into the controller as a validate method
    #if not file_location and not folder_location:
    #    print(f"Error: Neither file or folder locations were provided.")
    #    logger.error(f"Error: Neither file or folder locations were provided.")
    #    return

    #if file_location and folder_location:
    #    print(f"Error: Cannot have both file and folder locations.")
    #    logger.error(f"Error: Cannot have both file and folder locations.")
    #    return

    #if not mode:
    #    print(f"Error: No mode provided. --mode must be: contact_distance")
    #    logger.error(f"Error: No mode provided. --mode must be: contact_distance")
    #    return
    #else:
    #    if mode not in mode_list:
    #        print(f"Error: Mode {mode} is not supported. --mode must be: contact_distance")
    #        logger.error(f"Error: Mode {mode} is not supported. --mode must be: contact_distance")
    #        return
    #    else:
    #        print(f"Working in {mode} mode")
    #        logger.info(f"Working in {mode} mode")

    #if file_location:
    #    if not os.path.isfile(file_location):
    #        print(f"Error: The file {file_location} does not exist.")
    #        logger.error(f"Error: The file {file_location} does not exist.")
    #        return
    #    else:
    #        print(f"Processing input file {file_location}")
    #        logger.error(f"Processing input file {file_location}")

    #if folder_location:
    #    if not os.path.isdir(folder_location):
    #        print(f"Error: The folder {folder_location} does not exist.")
    #        logger.error(f"Error: The folder {folder_location} does not exist.")
    #        return
    #    else:
    #        print(f"Processing folder {folder_location}")
    #        logger.info(f"Processing folder {folder_location}")

    #if output_folder:
    #    if not os.path.isdir(output_folder):
    #        print(f"Error. The output folder {output_folder} does not exist. Create it first.")
    #        logger.error(f"Error. The output folder {output_folder} does not exist. Create it first.")
    #        return
    #else:
    #    print(f"Error: No output_folder specified. I don't know where to store results.")
    #    logger.error(f"Error: No output_folder specified. I don't know where to store results.")
    #    return

    #mode_decision(mode,
    #              ignore_num_start_res,
    #              ignore_num_end_res,
    #              output_folder,
    #              distance_cutoff,
    #              basis,
    #              theory,
    #              sequence_A,
    #              sequence_B,
    #              rotation_angle,
    #              output_file,
    #              renum_chains_list,
    #              file_location,
    #              folder_location)


#def mode_decision(mode: str,
#                  ignore_num_start_res: int,
#                  ignore_num_end_res: int,
#                  output_folder: str,
#                  distance_cutoff: float,
#                  basis: str,
#                  theory: str,
#                  sequence_A: str,
#                  sequence_B: str,
#                  rotation_angle: int,
#                  output_file: str,
#                  renum_chains_list: List[str] = None,
#                  file_location: str = None,
#                  folder_location: str = None) -> None:
#    logger = get_logger(__name__)
#    if mode == "contact_distance":
#        contact_controller = ContactController(file_location, folder_location, output_folder, distance_cutoff)
#        validated = contact_controller.validate_inputs(ignore_num_start_res, ignore_num_end_res, renum_chains_list)
#        if validated:
#            contact_controller.run_controller()
#        else:
#            print(f"Error: Unable to read inputs.")
#            logger.error(f"Error: Unable to read inputs.")
#    elif mode == "cap":
#        cap_controller = CapController(file_location, folder_location)
#        cap_controller.validate_inputs()
#        cap_controller.run_controller()
#    elif mode == "fmo":
#        fmo_controller = FMOController(basis, theory)
#        fmo_controller.validate_inputs(file_location, folder_location, output_folder)
#        fmo_controller.run_controller()
#    elif mode == "generate":
#        gen_helix_controller = GenHelixController(sequence_A, sequence_B, output_file)
#        validated = gen_helix_controller.validate_inputs()
#        if validated:
#            gen_helix_controller.run_controller()
#        else:
#            print(f"Error: Unable to read inputs for generating a helical dimer.")
#            logger.error(f"Error: Unable to read inputs for generating a helical dimer.")
#    elif mode == "rotation":
#        rotation_controller = RotationController(file_location, output_folder, rotation_angle)
#        rotation_controller.validate_inputs()
#        rotation_controller.run_controller()




if __name__ == '__main__':
    main()
