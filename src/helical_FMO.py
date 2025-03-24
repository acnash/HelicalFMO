# Main entry into HelicalFMO

import argparse
import os

from controllers.contact_controller import ContactController

def main() -> None:

    mode_list = ["contact_distance"]

    parser = argparse.ArgumentParser(description="Process TM domain PDB files or folder of files.")

    parser.add_argument("--file", type=str, help="Path to the file", required=False)
    parser.add_argument("--folder", type=str, help="Path to the folder", required=False)
    parser.add_argument("--mode", type=str, choices=["contact_distance"],
                        required=True, help="Select a mode: contact_distance.")
    parser.add_argument("--output_folder", type=str, required=True)
    parser.add_argument("--distance_cutoff", type=float, required=False, default=8)

    args = parser.parse_args()
    file_location = args.file
    folder_location = args.folder
    output_folder = args.output_folder
    distance_cutoff = args.distance_cutoff
    mode = args.mode

    if not file_location and not folder_location:
        print(f"Error: Neither file or folder locations were provided.")
        return

    if file_location and folder_location:
        print(f"Error: Cannot have both file and folder locations.")
        return

    if not mode:
        print(f"Error: No mode provided. --mode must be: contact_distance")
        return
    else:
        if mode not in mode_list:
            print(f"Error: Mode {mode} is not supported. --mode must be: contact_distance")
            return
        else:
            print(f"Working in {mode} mode")

    if file_location:
        if not os.path.isfile(file_location):
            print(f"Error: The file {file_location} does not exist.")
            return
        else:
            print(f"Processing input file {file_location}")

    if folder_location:
        if not os.path.isdir(folder_location):
            print(f"Error: The folder {folder_location} does not exist.")
            return
        else:
            print(f"Processing folder {folder_location}")

    if output_folder:
        if not os.path.isdir(output_folder):
            print(f"Error. The output folder {output_folder} does not exist. Create it first.")
            return
    else:
        print(f"Error: No output_folder specified. I don't know where to store results.")
        return

    mode_decision(mode, output_folder, distance_cutoff, file_location, folder_location)


def mode_decision(mode: str, output_folder: str, distance_cutoff: float, file_location: str = None, folder_location: str = None) -> None:
    if mode == "contact_distance":
        contact_controller = ContactController(file_location, folder_location, output_folder, distance_cutoff)
        validated = contact_controller.validate_inputs()
        if not validated:
            print(f"Error: Unable to read inputs.")
            return

        contact_controller.run_controller()




if __name__ == '__main__':
    main()
