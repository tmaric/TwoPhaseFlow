#!/bin/bash

grep -Po "Nnon\s*\K\d+" slurm-* | tee Nnon.dat > /dev/null

# Input and output file paths
input_file="Nnon.dat"
output_file="NnonPerStep.dat"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found."
   # exit 1
fi

if [ -e "$output_file" ]; then
    rm "$output_file"
   # exit 1
fi


# Read the input file line by line and perform the sum
echo "# Pressure equation solving times per time step" >> "$output_file"
echo -e "# times" >> "$output_file"

sum=0
count=0
while IFS= read -r line; do
    # Add the current number to the sum
    sum=$((sum + line))
    # Count the numbers read so far
    count=$((count + 1))

    # If we have read four numbers, calculate the sum and write to the output file
    if [ $((count % 4)) -eq 0 ]; then
        result=$((sum))
        step=$((count / 4))
        echo -e "$result" >> "$output_file"
        sum=0  # Reset the sum for the next four numbers
    fi
done < "$input_file"

rm "$input_file"

echo "Results written to '$output_file'."

