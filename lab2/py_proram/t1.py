import csv


class aohi:
    def __init__(self, row):
        self.data = row


def read_csv(file_path):
    data = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            data.append(row)
    return data

# Example usage
file_path = '../data/ewh.csv'
all_data = read_csv(file_path)


