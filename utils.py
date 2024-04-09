import numpy as np

def filter_array(arr):
    # Count occurrences of each number using numpy unique function
    unique_vals, counts = np.unique(arr, return_counts=True)

    # Create a dictionary to store counts of each unique value
    count_dict = dict(zip(unique_vals, counts))

    # Identify numbers that appear an even number of times
    even_occurrences = [num for num, count in count_dict.items() if count % 2 == 0]

    # Create a new filtered array
    filtered_arr = []
    seen_counts = {}

    # Iterate through the original array
    for num in arr:
        if num in even_occurrences:
            # Skip adding to filtered array if the number appears even times
            continue
        else:
            if num not in seen_counts:
                # Add the first occurrence of a number with odd times to filtered array
                filtered_arr.append(num)
                seen_counts[num] = 1
            # Otherwise, skip adding further occurrences of the same number

    return np.array(filtered_arr)


#variables from textbook
M1 = np.array([[0, 0, 5, 0, 0, 0, 0, 1], 
              [2, 0, 1, 0, 1, 1, 0, 1], 
              [0, 2, 0, 0, 0, 3, 0, 0], 
              [6, 2, 0, 0, 1, 0, 0, 0], 
              [1, 0, 0, 0, 0, 0, 0, 1], 
              [5, 0, 1, 0, 0, 2, 0, 0], 
              [0, 0, 2, 2, 0, 1, 0, 0]])
v1 = np.array([9398, 19095, 1964, 17078, 8077, 3397, 14262])
base = np.array([2, 3, 5, 7, 11, 13, 17, 19])