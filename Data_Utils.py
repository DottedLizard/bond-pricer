import csv
import re
import datetime


class DataUtils(object):
    """Class for data reading, writing, manipulation, and validation"""

    def __init__(self):
        pass

    @staticmethod
    def read_csv(file_name):

        csv_file = open(file_name, 'r', newline='\n')
        reader = csv.DictReader(csv_file)
        contents = []

        for row in reader:
            contents.append(row)

        csv_file.close()

        return contents

    @staticmethod
    def write_csv(file_name, data):

        csv_file = open(file_name, 'wt', newline='\n')
        writer = csv.writer(csv_file)

        if isinstance(data, dict):

            for key in data:
                DataUtils.print_row(writer, data[key])

        elif isinstance(data, list):

            for item in data:
                DataUtils.print_row(writer, item)

        else:
            DataUtils.print_row(writer, data)

        csv_file.close()

    @staticmethod
    def print_row(writer, item):

        if type(item) is str:
            writer.writerow([item])
        else:
            writer.writerow(item)

    @staticmethod
    def validate_in_list(item, lst, name):

        if item not in lst:
            raise ValueError(name + " must be one of " + str(lst))

    @staticmethod
    def validate_same_length(lst, name_lst):

        if not all([len(x) == len(y) for x, y in zip(lst, lst[1:])]):
            raise ValueError(str(name_lst) + " must be of same length")

    @staticmethod
    def validate_tenor(tenor):

        if re.match("^-?\d+[DMY]$", str(tenor)) is None:
            raise ValueError("Tenor must be in #D, #M, or #Y format")

    @staticmethod
    def validate_number(num, name):

        try:
            float(num)
        except ValueError:
            raise TypeError(name + " must be a number")

    @staticmethod
    def validate_date_or_tenor(date):

        if (type(date) is not datetime.date) and (re.match("^\d+[DMY]$", str(date)) is None):
            raise TypeError("Date must be of type datetime.date or a tenor string")

    @staticmethod
    def validate_type(inst, cls, name):

        if type(inst) is not cls:
            raise TypeError(name + " must be of type " + cls.__name__)
