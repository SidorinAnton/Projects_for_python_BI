# Created by Anton Sidorin (github - SidorinAnton),
# Anna Churkina (github - AnyaChurkina)
# and Margarita Komarova (github - Rita1612-GitHub)

import sys
import os
import time


def sequence_len(input_read: list, value: float) -> str:
    sequence = input_read[1]
    if len(sequence) >= value:
        return "pass"


def calculate_gc_content(input_read: list, value: list) -> str:
    sequence = input_read[1]
    lower_value = value[0]

    if len(value) == 2:
        upper_value = value[1]
        if lower_value > upper_value:
            lower_value, upper_value = upper_value, lower_value
    else:
        upper_value = None

    read_gc_bonds = ((sequence.upper().count("G") + sequence.upper().count("C")) / len(sequence)) * 100
    if upper_value is None:
        if read_gc_bonds >= lower_value:
            return "pass"
    else:
        if (read_gc_bonds >= lower_value) and (read_gc_bonds <= upper_value):
            return "pass"


def filter_the_reads(current_read: list, min_len: float, gc_content: list, failed=False) -> str:
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

    if sequence_len(current_read, min_len) == "pass" and calculate_gc_content(current_read, gc_content) == "pass":
        output_passed.write("\n".join(current_read))
        output_passed.write("\n")
        return "passed"
    else:
        if failed:
            output_failed.write("\n".join(current_read))
            output_failed.write("\n")
            return "failed"


if __name__ == "__main__":

    input_arguments = sys.argv
    arguments = dict()

    if len(input_arguments) == 1:
        print("Use '--help' to see the available options")
        exit()

    if "--help" in input_arguments:
        print(filter_the_reads.__doc__)
        exit()

    # Initialization of arguments

    # Another way of initialization
    # for val in input_arguments[1:]:
    #     # Works bad, if something is wrong with arguments (input errors)
    #     if val.startswith("--"):
    #         arguments[val] = []
    #         current_arg = val
    #     else:
    #         arguments[current_arg].append(val)

    if input_arguments[-1].endswith(".fastq"):
        file_with_reads = input_arguments.pop()
    else:
        raise FileNotFoundError("Input '.fastq' file is not defined. Check the type of file extension")

    if "--min_length" not in input_arguments:
        arguments["min_length"] = 0.0
        print("Parameter min_length - 0")
    else:
        arguments["min_length"] = float(input_arguments[input_arguments.index("--min_length") + 1])
        print(f"Parameter min_length - {arguments['min_length']}")

    if "--gc_bounds" not in input_arguments:
        arguments["gc_bounds"] = [0.0]
        print("Parameter gc_bounds - 0")
    else:
        arguments["gc_bounds"] = []
        for gc_val_index in range(input_arguments.index("--gc_bounds") + 1, len(input_arguments)):
            gc_value = input_arguments[gc_val_index]
            if gc_value.startswith("--"):
                break
            else:
                arguments["gc_bounds"].append(float(gc_value))
        print(f"Parameter gc_bounds - {arguments['gc_bounds']}")

    if "--output_base_name" in input_arguments:
        base_name = input_arguments[input_arguments.index("--output_base_name") + 1]
        print(f"Parameter output_base_name - {base_name}")
    else:
        base_name = file_with_reads[:len(file_with_reads) - 6]
        print(f"Parameter output_base_name - {base_name}")

    # If output files with the same base name exist
    if os.path.exists(f"{base_name}__passed.fastq"):
        print(f"{base_name}__passed.fastq exists")
        print("Do you want to rewrite this file? [y/n]")
        if input() == "y":
            os.remove(f"{base_name}__passed.fastq")
        else:
            print("Process finished")
            exit()
    if os.path.exists(f"{base_name}__failed.fastq"):
        print(f"{base_name}__failed.fastq exists")
        print("Do you want to rewrite this file? [y/n]")
        if input() == "y":
            os.remove(f"{base_name}__failed.fastq")
        else:
            print("Process finished")
            exit()

    # Main filtration
    start_time = time.time()
    with open(file_with_reads, "r") as input_file, open(f"{base_name}__passed.fastq", "a") as output_passed:
        read = []
        counter_passed = 0
        counter_failed = 0

        for line in input_file:
            if len(read) != 3:
                read.append(line.rstrip())
                continue
            else:
                read.append(line.rstrip())

                # Start of processing
                if "--keep_filtered" in input_arguments:
                    with open(f"{base_name}__failed.fastq", "a") as output_failed:
                        status = filter_the_reads(read, arguments["min_length"], arguments["gc_bounds"], failed=True)

                        if status == "passed":
                            counter_passed += 1
                        else:
                            counter_failed += 1
                else:
                    status = filter_the_reads(read, arguments["min_length"], arguments["gc_bounds"])

                    if status == "passed":
                        counter_passed += 1
                    else:
                        counter_failed += 1

            read = []

    if os.path.exists(f"{base_name}__failed.fastq"):
        print(f"Output files - {base_name}__passed.fastq and {base_name}__failed.fastq")
    else:
        print(f"Output file - {base_name}__passed.fastq")

    # End of processing and base report of filtration
    number_of_reads = counter_passed + counter_failed
    end_time = time.time()

    print()
    print("Running time", round(end_time - start_time, 3), "sec")
    print(f"{number_of_reads} reads filtered")
    print(f"{counter_passed} ({int(100 * counter_passed / number_of_reads)}%) pass filtration")
    print(f"{counter_failed} ({int(100 * counter_failed / number_of_reads)}%) fail filtration")
