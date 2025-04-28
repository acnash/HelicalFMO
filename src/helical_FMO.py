import argparse
import os
import yaml

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
        print(f"Making ../temp to store temporary files.")
        logger.info(f"Making ../temp to store temporary files.")
    else:
        print(f"../temp folder exists.")
        logger.info(f"../temp folder exists.")

    parser = argparse.ArgumentParser(description="Process TM domain PDB files or folder of files.")
    parser.add_argument("--input", type=str, required=True,
                        help=f"YAML input file. See the YAML files in ../examples folder for more information (and "
                             f"the GitHub README).")
    args = parser.parse_args()
    input_file = args.input

    with open(input_file, "r") as file:
        yaml_dict = yaml.safe_load(file)
        print(f"Opening YAML input file {input_file}.")
        logger.info(f"Opening YAML input file {input_file}.")

    for section in section_to_run:
        if section in yaml_dict:
            config_section = yaml_dict[section]

            if section == "generate":
                controller = GenHelixController()
            elif section == "rotation":
                controller = RotationController()
            elif section == "contact_distance":
                controller = ContactController()
            elif section == "cap":
                controller = CapController()
            elif section == "fmo":
                controller = FMOController()
            else:
                logger.warning(f"Warning: Operation {section} is not recognised.")
                print(f"Warning: Operation {section} is not recognised.")
                continue

            logger.info(f"Operation {section} is running.")
            print(f"Operation {section} is running.")

            valid = controller.validate_controller(config_section)
            if valid:
                print(f"Operation {section} is valid. Running controller...")
                logger.info(f"Operation {section} is valid. Running controller...")
                controller.run_controller()
            else:
                print(f"Error: Operation {section} is invalid. Aborting all operations.")
                logger.error(f"Operation {section} is invalid. Aborting all operations.")
                exit()


if __name__ == '__main__':
    main()
