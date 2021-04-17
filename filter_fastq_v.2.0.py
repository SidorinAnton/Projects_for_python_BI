"""
Filter_fastq v1.0.0
Usage:
    python3 filter_fastq.py [options] <input.fastq>
Available options:
--min_length <int> - minimal length of the sequence read to pass (Default = 0).
--keep_filtered - this flag should be set to save the not passed (failed the conditions) reads.
--gc_bounds <min> <max> - range of the GC bases percentage in a read. Use one value to specify only lower bound.
            Use two values to specify lower and upper bounds (Default = 0).
--output_base_name <str> - the basename of the output files (By default writes the base name of input
            without '.fastq')
Output files:
    output_base_name__passed.fastq
    output_base_name__failed.fastq (if 'keep_filtered' is defined)
"""
import os
import sys
import time
from typing import List


def check_arg_min_length(dict_args: dict, input_args: list,
                         invalid_argument: callable) -> dict:
    # Check --min_length
    if "--min_length" in input_args:
        index_of_value_after = input_args.index("--min_length") + 1
        min_len_value = input_args[index_of_value_after]

        if min_len_value.startswith("--"):
            raise invalid_argument()

        dict_args["--min_length"] = int(min_len_value)

    return dict_args


def check_arg_gc_content(dict_args: dict, input_args: list,
                         invalid_argument: callable) -> dict:
    # Check --gc_bounds
    if "--gc_bounds" in input_args:
        index_of_first_after_gc = input_args.index("--gc_bounds") + 1
        index_of_second_after_gc = input_args.index("--gc_bounds") + 2

        gc_first_value = input_args[index_of_first_after_gc]
        gc_second_value = input_args[index_of_second_after_gc]

        try:
            gc_first_value = float(gc_first_value)
        except ValueError:
            raise invalid_argument()
        else:
            dict_args["--gc_bounds"] = [gc_first_value]

        try:
            gc_second_value = float(gc_second_value)
        except ValueError:
            gc_second_value = None
        else:
            dict_args["--gc_bounds"].append(gc_second_value)

    return dict_args


def check_arg_keep_filtered(dict_args: dict, input_args: list) -> dict:
    if "--keep_filtered" in input_args:
        dict_args["--keep_filtered"] = True

    return dict_args


def check_arg_file_path(dict_args: dict, input_args: list,
                        invalid_argument: callable) -> dict:
    # Will work on Unix OS
    file_path = input_args[-1]
    file_name = file_path.split("/")[-1]

    if file_name.endswith(".fastq"):
        dict_args["file_path"] = file_path
        return dict_args

    raise invalid_argument()


def check_arg_output_base_name(dict_args: dict, input_args: list) -> dict:
    if "--output_base_name" in input_args:
        index_of_base_name = input_args.index("--output_base_name") + 1
        base_name = input_args[index_of_base_name]
    else:
        file_name = dict_args["file_path"].split("/")[-1]
        base_name, _ = file_name.split(".fastq")

    dict_args["--output_base_name"] = base_name
    return dict_args


def get_arguments(input_args: List[str]) -> dict:
    def invalid_argument():
        return TypeError(f"Your argument is incorrect.\nUse '--help' "
                         f"to see the available options")

    if len(input_args) == 1:
        print("Use '--help' to see the available options")
        exit()

    if "--help" in input_args:
        # print(filter_the_reads.__doc__)
        exit()

    dict_args = {
        "--min_length": 0,
        "--keep_filtered": False,
        "--gc_bounds": [0.0],
        "--output_base_name": None,
        "file_path": None
    }

    dict_args = check_arg_file_path(dict_args, input_args, invalid_argument)
    dict_args = check_arg_min_length(dict_args, input_args, invalid_argument)
    dict_args = check_arg_gc_content(dict_args, input_args, invalid_argument)
    dict_args = check_arg_keep_filtered(dict_args, input_args)
    dict_args = check_arg_output_base_name(dict_args, input_args)

    return dict_args


def sequence_len(input_read: list, value: int) -> bool:
    sequence = input_read[1]
    if len(sequence) >= value:
        return True
    return False


def calculate_gc_content(input_read: list, value: list) -> str:
    sequence = input_read[1]
    lower_value = value[0]

    if len(value) == 2:
        upper_value = value[1]
        if lower_value > upper_value:
            lower_value, upper_value = upper_value, lower_value
    else:
        upper_value = None

    read_gc_bonds = ((sequence.upper().count("G") + sequence.upper().count("C"))
                     / len(sequence)) * 100
    if upper_value is None:
        return read_gc_bonds >= lower_value
    else:
        return (read_gc_bonds >= lower_value) and (read_gc_bonds <= upper_value)


def check_file_exists(output_base_name):
    if os.path.exists(f"{output_base_name}__passed.fastq"):
        print(f"{output_base_name}__passed.fastq exists")
        print("Do you want to rewrite this file? [y/n]")
        if input() == "y":
            os.remove(f"{output_base_name}__passed.fastq")
        else:
            print("Process finished")
            exit()
    if os.path.exists(f"{output_base_name}__failed.fastq"):
        print(f"{output_base_name}__failed.fastq exists")
        print("Do you want to rewrite this file? [y/n]")
        if input() == "y":
            os.remove(f"{output_base_name}__failed.fastq")
        else:
            print("Process finished")
            exit()


def write_processed_files(dict_args):
    """

    :param dict_:
    :return:
    """

    min_len = dict_args["--min_length"]
    gc = dict_args["--gc_bounds"]
    file_path = dict_args["file_path"]
    keep_filtered = dict_args["--keep_filtered"]
    output_base_name = dict_args["--output_base_name"]
    check_file_exists(output_base_name)
    output_file_passed = output_base_name + "__passed.fastq"

    with open(file_path, "r") as fastq_file:
        read = []
        counter_passed = 0
        counter_failed = 0

        for line in fastq_file:
            if len(read) != 3:
                read.append(line.rstrip())
                continue
            else:
                read.append(line.rstrip())

            len_status = sequence_len(read, min_len)
            gc_status = calculate_gc_content(read, gc)

            read_status = len_status and gc_status
            with open(output_file_passed, "a") as output_passed:
                if read_status:
                    output_passed.write("\n".join(read))
                    output_passed.write("\n")
                    counter_passed += 1
                else:
                    counter_failed += 1
                    if keep_filtered:
                        output_file_failed = output_base_name + "__failed.fastq"
                        with open(output_file_failed, "a") as output_failed:
                            output_failed.write("\n".join(read))
                            output_failed.write("\n")
            read = []
    return counter_passed, counter_failed


if __name__ == "__main__":
    start_time = time.time()
    input_arguments = sys.argv
    dict_arguments = get_arguments(input_arguments)
    counter_passed, counter_failed = write_processed_files(dict_arguments)
    number_of_reads = counter_passed + counter_failed
    end_time = time.time()

    print()
    print("Running time", round(end_time - start_time, 3), "sec")
    print(f"{number_of_reads} reads filtered")
    print(f"{counter_passed} ({int(100 * counter_passed / number_of_reads)}%) "
          f"pass filtration")
    print(f"{counter_failed} ({int(100 * counter_failed / number_of_reads)}%) "
          f"fail filtration")
