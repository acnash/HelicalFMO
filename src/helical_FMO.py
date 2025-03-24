#Main entry into HElicalFMO

import argparse
import os

def main():

    parser = argparse.ArgumentParser(description="Process TM domain PDB files or folder of files.")

    parser.add_argument("--file", type=str, help="Path to the file", required=False)
    parser.add_argument("--folder", type=str, help="Path to the folder", required=False)

    args = parser.parse_args()
    file_location = args.file
    folder_location = args.folder

    if not file_location and not folder_location:
        print(f"Error: Neither file or folder locations were provided.")
        return

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



if __name__ == "__main__":
    main()