import unittest
from Fastq_filtrator import *
import os


class Test_check_arg_min_length(unittest.TestCase):

    def test_correct_number(self):
        dict_result = check_arg_min_length({}, ["--min_length", '2'],
                                           ValueError)
        actual_result = dict_result["--min_length"]
        self.assertEqual(actual_result, 2)

    def test_absence_of_value_of_argument(self):
        self.assertRaises(RecursionError, check_arg_min_length, {},
                          ["--min_length", '--gc_bounds'], RecursionError)

    def test_incorrect_name_of_argument(self):
        dict_result = check_arg_min_length({}, ["--min_len", '2'], ValueError)
        actual_result = len(dict_result)
        self.assertEqual(actual_result, 0)


class Test_check_arg_gc_content(unittest.TestCase):

    def test_correct_first_number(self):
        dict_result = check_arg_gc_content({},
                                           ["--gc_bounds", '2', 'something'],
                                           ValueError)
        actual_result = dict_result["--gc_bounds"]
        self.assertEqual(actual_result, [2])

    def test_correct_second_number(self):
        dict_result = check_arg_gc_content({}, ["--gc_bounds", '2', '3'],
                                           ValueError)
        actual_result = dict_result["--gc_bounds"]
        self.assertEqual(actual_result, [2, 3])

    def test_absence_of_value_of_argument(self):
        self.assertRaises(RecursionError, check_arg_gc_content, {},
                          ["--gc_bounds", '', '--min_length'], RecursionError)

    def test_incorrect_name_of_argument(self):
        dict_result = check_arg_gc_content({}, ["--gc", '2', '3'], ValueError)
        actual_result = len(dict_result)
        self.assertEqual(actual_result, 0)


class Test_check_arg_keep_filtered(unittest.TestCase):

    def test_absence_of_assign_of_argument(self):
        dict_result = check_arg_keep_filtered({}, ["--keep_filtered"])
        actual_result = dict_result["--keep_filtered"]
        self.assertTrue(actual_result)

    def test_incorrect_name_of_argument(self):
        dict_result = check_arg_keep_filtered({}, ["--keep_filt"])
        actual_result = len(dict_result)
        self.assertEqual(actual_result, 0)


class Test_check_arg_file_path(unittest.TestCase):

    def test_type_of_file(self):
        self.assertRaises(RecursionError, check_arg_file_path, {},
                          ['/python_HW/file'], RecursionError)

    def test_correct_filepath(self):
        dict_result = check_arg_file_path({}, ["/python_HW/file.fastq"],
                                          ValueError)
        actual_result = dict_result['file_path']
        self.assertEqual(actual_result, '/python_HW/file.fastq')


class Test_check_arg_output_base_name(unittest.TestCase):

    def test_absence_base_name(self):
        dict_result = check_arg_output_base_name({"file_path": "file.fastq"},
                                                 [])
        expected_result = "file"
        actual_result = dict_result["--output_base_name"]
        self.assertEqual(actual_result, expected_result)

    def test_presence_base_name(self):
        dict_result = check_arg_output_base_name({"file_path": "file.fastq"},
                                                 ["--output_base_name",
                                                  'fastq_name'])
        expected_result = "fastq_name"
        actual_result = dict_result["--output_base_name"]
        self.assertEqual(actual_result, expected_result)

    def test_incorrect_name_of_argument(self):
        dict_result = check_arg_output_base_name({"file_path": "file.fastq"},
                                                 ["--output_name",
                                                  'fastq_name'])
        expected_result = "file"
        actual_result = dict_result["--output_base_name"]
        self.assertEqual(actual_result, expected_result)


class Test_get_arguments(unittest.TestCase):

    def test_input_arguments(self):
        input_args = [
            ["filtrator.py", "--min_length", "2", "--gc_bounds", "3",
             "file.fastq"],
            ["filtrator.py", "--gc_bounds", "3", "20", "file.fastq"],
            ["filtrator.py", "--keep_filtered", "/path/to/file.fastq"],
            ["filtrator.py", "/path/to/file.fastq"],
            ["filtrator.py", "--output_base_name", "name",
             "/path/to/file.fastq"]
        ]

        result = [
            {"--min_length": 2, "--gc_bounds": [3], "file_path": "file.fastq",
             "--output_base_name": "file",
             "--keep_filtered": False},

            {"--min_length": 0, "--gc_bounds": [3, 20],
             "file_path": "file.fastq", "--output_base_name": "file",
             "--keep_filtered": False},

            {"--min_length": 0, "--gc_bounds": [0],
             "file_path": "/path/to/file.fastq", "--output_base_name": "file",
             "--keep_filtered": True},

            {"--min_length": 0, "--gc_bounds": [0],
             "file_path": "/path/to/file.fastq", "--output_base_name": "file",
             "--keep_filtered": False},

            {"--min_length": 0, "--gc_bounds": [0],
             "file_path": "/path/to/file.fastq", "--output_base_name": "name",
             "--keep_filtered": False},

        ]

        for list_args, expected_res in zip(input_args, result):
            actual_result = get_arguments(list_args)
            with self.subTest(i=(list_args, expected_res)):
                self.assertDictEqual(actual_result, expected_res)


class Test_sequence_len(unittest.TestCase):
    def test_pass_length_of_seq(self):
        actual_result = sequence_len(['', 'ATGC', '', ''], 2)
        self.assertTrue(actual_result)

    def test_failed_length_of_seq(self):
        actual_result = sequence_len(['', 'ATGC', '', ''], 5)
        self.assertFalse(actual_result)


class Test_calculate_gc_content(unittest.TestCase):
    def test_pass_gc_with_one_value(self):
        actual_result = calculate_gc_content(['', 'ATGC', '', ''], [20])
        self.assertTrue(actual_result)

    def test_failed_gc_with_one_value(self):
        actual_result = calculate_gc_content(['', 'ATGC', '', ''], [80])
        self.assertFalse(actual_result)

    def test_pass_gc_with_two_value(self):
        actual_result = calculate_gc_content(['', 'ATGC', '', ''], [20, 60])
        self.assertTrue(actual_result)

    def test_failed_gc_with_two_value(self):
        actual_result = calculate_gc_content(['', 'ATGC', '', ''], [70, 90])
        self.assertFalse(actual_result)


class Test_write_processed_files(unittest.TestCase):
    def setUp(self) -> None:
        with open("name.fastq", "w") as file:
            a = [
                "", "ATGCATGC", "", "",  # Len = 8, GC = 50
                "", "ATGC", "", "",  # Len = 4, GC = 50
                "", "GGGG", "", "",  # Len = 4, GC = 100
                "", "AAAA", "", ""  # Len = 4, GC = 0
            ]

            for i in a:
                file.write(i + "\n")

    def tearDown(self) -> None:
        os.remove("name.fastq")  # Input

        #     os.remove("name__failed.fastq")  # No output_base_name
        #     os.remove("test_bn__failed.fastq")  # output_base_name = "test_bn"

    def test_correct_passed_with_basename(self):
        input_arguments = {"--min_length": 2, "--gc_bounds": [0],
                           "file_path": "name.fastq",
                           "--output_base_name": "test_bn",
                           "--keep_filtered": False}
        write_processed_files(input_arguments)

        expected_result_passed = [
            '', 'ATGCATGC', '', '',
            '', 'ATGC', '', '',
            '', 'GGGG', '', '',
            '', 'AAAA', '', ''
        ]

        with open("test_bn__passed.fastq") as test:
            actual_result = []
            for line in test:
                actual_result.append(line.strip())

        os.remove("test_bn__passed.fastq")
        self.assertListEqual(expected_result_passed, actual_result)

    def test_correct_passed_without_basename(self):
        input_arguments = {"--min_length": 2, "--gc_bounds": [0],
                           "file_path": "name.fastq",
                           "--output_base_name": "name",
                           "--keep_filtered": False}
        write_processed_files(input_arguments)

        expected_result_passed = [
            '', 'ATGCATGC', '', '',
            '', 'ATGC', '', '',
            '', 'GGGG', '', '',
            '', 'AAAA', '', ''
        ]

        with open("name__passed.fastq") as test:
            actual_result = []
            for line in test:
                actual_result.append(line.strip())

        os.remove("name__passed.fastq")
        self.assertListEqual(expected_result_passed, actual_result)

    def test_correct_failed_with_basename(self):
        input_arguments = {"--min_length": 100, "--gc_bounds": [0],
                           "file_path": "name.fastq",
                           "--output_base_name": "test_bn",
                           "--keep_filtered": False}
        write_processed_files(input_arguments)

        expected_result_failed = []

        with open("test_bn__passed.fastq") as test:
            actual_result = []
            for line in test:
                actual_result.append(line.strip())

        os.remove("test_bn__passed.fastq")
        self.assertListEqual(expected_result_failed, actual_result)

    def test_correct_failed_without_basename(self):
        input_arguments = {"--min_length": 100, "--gc_bounds": [0],
                           "file_path": "name.fastq",
                           "--output_base_name": "name",
                           "--keep_filtered": False}
        write_processed_files(input_arguments)

        expected_result_failed = []

        with open("name__passed.fastq") as test:
            actual_result = []
            for line in test:
                actual_result.append(line.strip())

        os.remove("name__passed.fastq")
        self.assertListEqual(expected_result_failed, actual_result)

    def test_correct_passed_and_failed_with_basename(self):
        input_arguments = {"--min_length": 5, "--gc_bounds": [0],
                           "file_path": "name.fastq",
                           "--output_base_name": "test_bn",
                           "--keep_filtered": True}
        write_processed_files(input_arguments)

        expected_result_passed = [
            '', 'ATGCATGC', '', '']

        expected_result_failed = [
            '', 'ATGC', '', '',
            '', 'GGGG', '', '',
            '', 'AAAA', '', '']

        with open("test_bn__passed.fastq") as test:
            actual_result_passed = []
            for line in test:
                actual_result_passed.append(line.strip())

        with open("test_bn__failed.fastq") as test:
            actual_result_failed = []
            for line in test:
                actual_result_failed.append(line.strip())

        os.remove("test_bn__passed.fastq")
        os.remove("test_bn__failed.fastq")
        self.assertListEqual(expected_result_passed, actual_result_passed)
        self.assertListEqual(expected_result_failed, actual_result_failed)

    def test_correct_passed_and_failed_without_basename(self):
        input_arguments = {"--min_length": 5, "--gc_bounds": [0],
                           "file_path": "name.fastq",
                           "--output_base_name": "name",
                           "--keep_filtered": True}
        write_processed_files(input_arguments)

        expected_result_passed = [
            '', 'ATGCATGC', '', '']

        expected_result_failed = [
            '', 'ATGC', '', '',
            '', 'GGGG', '', '',
            '', 'AAAA', '', '']

        with open("name__passed.fastq") as test:
            actual_result_passed = []
            for line in test:
                actual_result_passed.append(line.strip())

        with open("name__failed.fastq") as test:
            actual_result_failed = []
            for line in test:
                actual_result_failed.append(line.strip())

        os.remove("name__passed.fastq")
        os.remove("name__failed.fastq")
        self.assertListEqual(expected_result_passed, actual_result_passed)
        self.assertListEqual(expected_result_failed, actual_result_failed)


if __name__ == '__main__':
    unittest.main()
