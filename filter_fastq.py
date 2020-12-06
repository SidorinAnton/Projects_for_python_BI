
import sys
import os.path
import time


def sequence_len(input_read: list, value):
    # print("In function length_of_the_sequence -", read)
    sequence = input_read[1]
    if len(sequence) >= float(value):
        return "pass"


def calculate_gc_content(input_read: list, value: list):
    sequence = input_read[1]
    lower_value = float(value[0])

    if len(value) == 2:
        upper_value = float(value[1])
    else:
        upper_value = None

    read_gc_bonds = ((sequence.upper().count("G") + sequence.upper().count("C")) / len(sequence)) * 100
    if upper_value is None:
        if read_gc_bonds >= lower_value:
            return "pass"
    else:
        if (read_gc_bonds >= lower_value) and (read_gc_bonds <= upper_value):
            return "pass"


def filter_the_reads(current_read: list, min_len, gc_content: list, failed=False):
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

    start_time = time.time()

    input_arguments = sys.argv
    file_with_reads = input_arguments.pop()
    arguments = dict()

    for val in input_arguments[1:]:
        # Works bad, if something is wrong with arguments (input errors)
        if val.startswith("--"):
            arguments[val] = []
            current_arg = val
        else:
            arguments[current_arg].append(val)

    if "--min_length" not in arguments:
        raise NameError("Argument '--min_length' should be defined")

    if "--gc_bounds" not in arguments:
        arguments["--gc_bounds"] = [0]  # GC content will be greater, than zero

    if "--output_base_name" in arguments:
        base_name = arguments["--output_base_name"][0]
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
                print("Processing read", read[0])

                if "--keep_filtered" in arguments:
                    with open(f"{base_name}__failed.fastq", "a") as output_failed:
                        status = filter_the_reads(read, arguments["--min_length"][0], arguments["--gc_bounds"], failed=True)

                        if status == "passed":
                            counter_passed += 1
                        else:
                            counter_failed += 1
                else:
                    status = filter_the_reads(read, arguments["--min_length"][0], arguments["--gc_bounds"])

                    if status == "passed":
                        counter_passed += 1
                    else:
                        counter_failed += 1

            read = []
    number_of_reads = counter_passed + counter_failed
    end_time = time.time()

    print()
    print("Running time", round(end_time - start_time, 3))
    print(f"{number_of_reads} reads filtered")
    print(f"{counter_passed} ({int(100 * counter_passed / number_of_reads)}%) pass filtration")
    print(f"{counter_failed} ({int(100 * counter_failed / number_of_reads)}%) fail filtration")

