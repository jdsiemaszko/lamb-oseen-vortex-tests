import numpy as np
import csv
import json

arr = np.loadtxt("initial_data.csv", delimiter = ' ')



# Define the input and output file paths
json_file = 'particles.json'

# Initialize an empty list to store particle data
particle_data = []

for row in arr:
    # Define a particle object for each row
    particle = {
        "center": [float(row[0]), float(row[1])],
        "enabled": True,
        "strength": float(row[2]),
        "type": "single particle"
    }
    # Append the particle object to the list
    particle_data.append(particle)

# Create a dictionary with the "flowstructures" key and the particle data
particle_json = {"flowstructures": particle_data}

# Write the JSON data to the output file
with open(json_file, 'w') as json_file:
    json.dump(particle_json, json_file, indent=2)

print(f"JSON file '{json_file}' has been created with {len(particle_data)} particles.")