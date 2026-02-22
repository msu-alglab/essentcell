import argparse
import itertools
import math
import os

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Arguments for doublet filtering program.")
    parser.add_argument('inputFile', type=str, help="Input file to filter")
    parser.add_argument('cells_to_remove', type=int, default=1,
                        help="Number of Cells that needs to be removed after filtering")
    parser.add_argument('-output_folder', type=str, default="filtered",
                        help="The output folder name for filtered files")
    args = parser.parse_args()
    print(f"Input file: {args.inputFile}")
    print(f"Cells to remove: {args.cells_to_remove}")
    print(f"Output folder: {args.output_folder}")
    print("**********************************")

    input_file_name = args.inputFile
    output_folder = args.output_folder
    cells_to_remove = args.cells_to_remove

    df = pd.read_csv(args.inputFile)
    column_list = df.columns.tolist()

    df.head()
    D = df.to_numpy()
    rows, cols = D.shape
    rank_dict = compute_rank(D)

    print(f"Rank dictionary: {rank_dict}")
    sorted__pairs_of_ranks = sorted(rank_dict.items(), key=lambda x: x[1], reverse=True)
    if cells_to_remove > rows:
        cells_to_remove = rows
    number_of_rows_to_delete = cells_to_remove
    print(f"Number of rows to delete: {number_of_rows_to_delete}")

    indices_to_remove = [i for i, j in sorted__pairs_of_ranks[:number_of_rows_to_delete]]
    print(f"indices_to_remove: {indices_to_remove}")

    new_D = np.delete(D, indices_to_remove, axis=0)
    print(f"new_D.shape: {new_D.shape}")

    output_file_name = input_file_name[:-4].split('/')[-1]  # Extract the file name without the extension
    filtered_df = pd.DataFrame(new_D, columns=column_list)

    if not os.path.exists(f'{output_folder}/{output_file_name}'):
        os.makedirs(f'{output_folder}/{output_file_name}')
    filtered_df.to_csv(f"{output_folder}/{output_file_name}/{output_file_name}_{number_of_rows_to_delete}_filtered.csv",
                       index=False)


def compute_rank(D) -> dict[int, int]:
    rows, cols = D.shape
    rank_dict = {i: 0 for i in range(rows)}
    for row in range(rows):
        rank_dict[row] = compute_rank_for_k(D, row)
    print(rank_dict)
    return rank_dict


def compute_rank_for_k(data: np.ndarray, k: int) -> int:
    rows, cols = data.shape
    numbers = [i for i in range(rows) if i != k]
    # print(numbers)
    rank = 0
    for i, j in itertools.combinations(numbers, 2):
        rank += compute_doublety_value(data[i], data[j], data[k])
    return rank


def compute_doublety_value(i_array, j_array, k_array) -> int:
    assert i_array.shape[0] == j_array.shape[0]
    assert i_array.shape[0] == k_array.shape[0]
    cols = i_array.shape[0]
    count = 0;
    numbers = [i for i in range(0, cols)]
    for p, q in itertools.combinations(numbers, 2):
        if (i_array[p] == j_array[q] == 0) and (i_array[q] == j_array[p] == k_array[p] == k_array[q] == 1):
            count += 1
        elif (i_array[q] == j_array[p] == 0) and (i_array[p] == j_array[q] == k_array[p] == k_array[q] == 1):
            count += 1
    return count


def test_compute_rank():
    assert compute_rank_for_k(np.array([[0, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 1]]), 1) == 0
    assert compute_rank_for_k(np.array([[0, 0, 0, 1], [0, 1, 0, 0], [0, 1, 0, 1]]), 2) == 1


def test_full_compute_rank():
    D = np.array([[0, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 1]])
    assert compute_rank_for_k(D, 0) == 0
    assert compute_rank_for_k(D, 1) == 0
    assert compute_rank_for_k(D, 2) == 1


if __name__ == '__main__':
    main()
