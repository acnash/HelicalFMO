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

    parser = argparse.ArgumentParser(description="Process TM domain PDB files or folder of files.")
    parser.add_argument("--input", type=str, required=True, help="YAML input file.")
    args = parser.parse_args()
    input_file = args.input

    with open(input_file, "r") as file:
        yaml_dict = yaml.safe_load(file)

    for section in section_to_run:
        logger.info(f"Checking operation {section}...")
        print(f"Checking operation {section}...")
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

            controller.validate_controller(config_section)
            controller.run_controller()


if __name__ == '__main__':
    main()
