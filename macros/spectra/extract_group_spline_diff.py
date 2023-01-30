# Usage: python extract_group_spline_diff.py OutputFile.txt (--mode_option)
# A simple python script that takes a txt output file from 'plot_sigma_var_all_groups.py" and returns all params (+modes) that show difference to command line.
# Getting an output.txt file from this can be done by running: Script | tee Output.txt
# Any complaints, send to Jaiden
# Jan 2023
#
import argparse

def extract_string(sentence, mode_option):
    # Define the start of the "mode" and "param" substrings
    start_mode = "mode = "
    start_param = "For param = "
    end = " ; "

    # Find the start index of the "mode" substring
    start_index_mode = sentence.find(start_mode)
    # If not found, return None
    if start_index_mode == -1:
        return None
    # Add the length of "mode = " to the start index to find the start of the "mode" value
    start_index_mode += len(start_mode)
    # Find the end index of the "mode" substring
    end_index_mode = sentence.find(end, start_index_mode)
    # Extract the "mode" value
    mode = sentence[start_index_mode:end_index_mode]

    # Find the start index of the "param" substring
    start_index_param = sentence.find(start_param)
    # If not found, return None
    if start_index_param == -1:
        return None
    # Add the length of "For param = " to the start index to find the start of the "param" value
    start_index_param += len(start_param)
    # Find the end index of the "param" substring
    end_index_param = sentence.find(end, start_index_param)
    # Extract the "param" value
    param = sentence[start_index_param:end_index_param]

    # Return either the "param" value only, or the combined "param" and "mode" values
    return param if not mode_option else mode + " " + param

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Input file name')
    parser.add_argument('--mode_option', action='store_true', help='Boolean flag to return only "param" if False')
    args = parser.parse_args()
    filename = args.filename
    mode_option = args.mode_option
    strings_seen = set()
    # Open the input file
    with open(filename) as f:
        # Read all the lines in the file
        lines = f.readlines()
        for line in lines:
            result = extract_string(line, mode_option)
            # If a string was extracted and it hasn't been seen before
            if result and result not in strings_seen:
                strings_seen.add(result)
                print(result)

