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

    args = parser.parse_args()
    file_location = args.file
    folder_location = args.folder
    mode = args.mode

    if not file_location and not folder_location:
        print(f"Error: Neither file or folder locations were provided.")
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

    mode_decision(mode)

def mode_decision(mode: str) -> None:
    if mode == "contact_distance":
        contact_controller = ContactController()




if __name__ == '__main__':
    main()
