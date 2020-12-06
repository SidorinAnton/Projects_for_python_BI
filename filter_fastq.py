# Created by Anton Sidorin (github - SidorinAnton),
# Anna Churkina (github - AnyaChurkina), Margarita Komarova (github - Rita1612-GitHub)


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
    :param current_read:
    :param min_len:
    :param gc_content:
    :param failed:
    :return:
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

    print("Use '--help' to see ....")
    input_arguments = sys.argv

    if "--help" in input_arguments:
        print(filter_the_reads.__doc__)

    start_time = time.time()

    file_with_reads = input_arguments.pop()
    arguments = dict()

    # for val in input_arguments[1:]:
    #     # Works bad, if something is wrong with arguments (input errors)
    #     if val.startswith("--"):
    #         arguments[val] = []
    #         current_arg = val
    #     else:
    #         arguments[current_arg].append(val)

    if "--min_length" not in input_arguments:
        raise NameError("Argument '--min_length' should be defined")
    else:
        value_of_min_len = input_arguments[input_arguments.index("--min_length") + 1]
        arguments["min_length"] = float(value_of_min_len)




    if "--gc_bounds" not in input_arguments:
        arguments["gc_bounds"] = [0]  # GC content will be greater, than zero
    else:
        arguments["gc_bounds"] = []
        index_of_gc_value = input_arguments.index("--gc_bounds")

        for gc_val_index in range(index_of_gc_value + 1, len(input_arguments)):
            gc_value = input_arguments[gc_val_index]
            if gc_value.startswith("--"):
                break
            else:
                arguments["gc_bounds"].append(float(gc_value))



    if "--output_base_name" in input_arguments:
        index_of_base_name = input_arguments.index("--output_base_name") + 1
        base_name = input_arguments[index_of_base_name]
    else:
        base_name = file_with_reads[:len(file_with_reads) - 6]



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



    with open(file_with_reads, "r") as input_file, open(f"{base_name}__passed.fastq", "a") as output_passed:
        read = []
        counter_passed = 0
        counter_failed = 0

        for smt in input_file:
            if len(read) != 3:
                read.append(smt.rstrip())
                continue
            else:
                read.append(smt.rstrip())
                # print("Processing read", read[0])

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
    number_of_reads = counter_passed + counter_failed
    end_time = time.time()

    print()
    print("Running time", round(end_time - start_time, 3), "sec")
    print(f"{number_of_reads} reads filtered")
    print(f"{counter_passed} ({int(100 * counter_passed / number_of_reads)}%) pass filtration")
    print(f"{counter_failed} ({int(100 * counter_failed / number_of_reads)}%) fail filtration")

