import matplotlib.pyplot as plt

# Function to read a and p points from the file
def read_points_from_file(file_path):
    a_points = []  # to store a points
    p_points = []  # to store p points
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith('a'):
                # Extracting the coordinates and converting them to float
                parts = line.split()
                a_points.append((float(parts[1]), float(parts[2])))
            elif line.startswith('p'):
                parts = line.split()
                p_points.append((float(parts[1]), float(parts[2])))

    return a_points, p_points

# Function to plot the points
def plot_and_save(a_points, p_points, output_path):
    # Extract x and y coordinates
    a_x, a_y = zip(*a_points) if a_points else ([], [])
    p_x, p_y = zip(*p_points) if p_points else ([], [])

    # Create plot
    plt.figure(figsize=(6, 6))  # Square plot
    plt.scatter(a_x, a_y, color='red', label='Moveable Cells', s=5)
    plt.scatter(p_x, p_y, color='blue', label='I/O Pads', s=5)

    # Adding title and labels with a specified fontsize
    plt.title("Placement before spreading", fontsize=20)
    plt.xlabel("X", fontsize=14)
    plt.ylabel("Y", fontsize=14)

    # Position the legend outside of the plot area on the bottom
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2)

    # Adjust the plot area to include the legend and labels
    plt.tight_layout()

    # Save the plot as a file
    plt.savefig(output_path, dpi=400)
    plt.close()  # Close the figure to prevent display if running in a non-interactive environment

# Replace with your actual file paths
input_file_path = 'output/coordinates.txt' 
output_file_path = 'output/plot.jpg'

# Now we can call the functions with the file paths
a_points, p_points = read_points_from_file(input_file_path)
plot_and_save(a_points, p_points, output_file_path)
