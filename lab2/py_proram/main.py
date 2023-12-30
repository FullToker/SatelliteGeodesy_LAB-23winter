### by Yushuo, 29.12.2023
### satellite geodesy lab2

import scipy.io
import pandas as pd
import csv

def mat_to_csv(mat_file_path, csv_file_path, name):
    # Load the MATLAB .mat file
    mat_data = scipy.io.loadmat(mat_file_path)

    # Extract the data from the MATLAB file (modify this based on your file structure)
    data_to_export = mat_data[name]

    # Write the data to a CSV file
    with open(csv_file_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)

        # Write header if needed
        # csvwriter.writerow(['Column1', 'Column2', ...])

        # Write data rows
        for rows in data_to_export:
                for row in rows:
                    csvwriter.writerow(row)

def mat_to_csv2(mat_file, csv_file):
    # Load MATLAB file
    mat_data = scipy.io.loadmat(mat_file)

    # Convert to DataFrame
    df = pd.DataFrame(mat_data)

    # Save DataFrame to CSV
    df.to_csv(csv_file, index=False)

if __name__ == "__main__":
    # Replace 'input.mat' with the path to your MATLAB file
    input_mat_file = './data/lab2_data.mat'

    # Replace 'output.csv' with the desired CSV file name
    output_csv_file = './data/lab2.csv'

    def1='./data/def.csv'
    ewh='./data/ewh.csv'
    mat_to_csv(input_mat_file, def1 , 'aohi_def')
    mat_to_csv(input_mat_file, ewh , 'aohi_ewh')


    # mat_to_csv2(input_mat_file, output_csv_file)
