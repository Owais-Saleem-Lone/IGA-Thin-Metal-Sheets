import numpy as np

# Your original array
original_array = np.array([0, 5, 10, 15, 20])

# Number of values to insert between each interval
num_values_between = 3

# Initialize an empty list to store the new values
new_values = []

# Iterate through each interval in the original array
for i in range(len(original_array) - 1):
    # Generate 'num_values_between' values between the current and next interval
    values_between = np.linspace(original_array[i], original_array[i + 1], num_values_between + 2)[1:-1]
    
    # Append the generated values to the new_values list
    # new_values.extend(values_between)
    new_values.extend([original_array[i]] + list(values_between) + [original_array[i + 1]])

# Convert the list to a NumPy array
new_array = np.unique(np.array(new_values))

# Print the result
print(new_array)
