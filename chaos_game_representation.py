import matplotlib.pyplot as plt
import numpy as np

def plot_cgr(sequences, legend_labels=None):
    # Mapping nucleotides to coordinates
    mapping = {
        'A': (0, 0),
        'C': (1, 0),
        'G': (1, 1),
        'T': (0, 1)
    }

    # Generate a unique color for each sequence
    num_sequences = len(sequences)
    colors = plt.cm.tab10(np.linspace(0, 1, num_sequences))

    # Set up the plot
    plt.figure(figsize=(8, 8))
    plt.title("Chaos Game Representation")
    plt.xlabel("X")
    plt.ylabel("Y")

    # Iterate over each sequence
    for i, sequence in enumerate(sequences):
        # Initial position
        current_pos = (0.5, 0.5)

        # List to store the X and Y coordinates
        x_coords = [current_pos[0]]
        y_coords = [current_pos[1]]

        # Generate CGR for the current sequence
        for j, nucleotide in enumerate(sequence):
            # Get the corresponding coordinates for the nucleotide
            new_pos = mapping.get(nucleotide)

            if new_pos:
                # Calculate the midpoint between the current position and the new position
                current_pos = (
                    (current_pos[0] + new_pos[0]) / 2,
                    (current_pos[1] + new_pos[1]) / 2
                )

                # Append the X and Y coordinates
                x_coords.append(current_pos[0])
                y_coords.append(current_pos[1])

                # Add number label to each dot
                plt.text(current_pos[0], current_pos[1], str(j+1), color=colors[i], fontsize=8)

        # Add scatter plot for the current sequence with the corresponding color
        label = legend_labels[i] if legend_labels else f"Sequence {i+1}"
        plt.scatter(x_coords, y_coords, label=label, color=colors[i], s=5)

        # Connect the dots of the sequence
        plt.plot(x_coords, y_coords, color=colors[i], linewidth=1, alpha=0.6)

    # Draw the dividing lines and label each quadrant
    plt.axvline(x=0.5, linestyle='dashed', color='grey')
    plt.axhline(y=0.5, linestyle='dashed', color='grey')

    # Add the letters 'A', 'T', 'C', 'G' inside their respective quadrants
    plt.text(0.02, 0.02, 'A', fontsize=12, color="grey" , fontweight='bold', transform=plt.gca().transAxes)
    plt.text(0.98, 0.02, 'T', fontsize=12, color="grey" , fontweight='bold', ha='right', transform=plt.gca().transAxes)
    plt.text(0.02, 0.98, 'C', fontsize=12, color="grey" , fontweight='bold', va='top', transform=plt.gca().transAxes)
    plt.text(0.98, 0.98, 'G', fontsize=12, color="grey" , fontweight='bold', ha='right', va='top', transform=plt.gca().transAxes)

    # Change the background color of the upper right and lower left quadrants to light gray
    plt.gca().add_patch(plt.Rectangle((0, 0), 0.5, 0.5, facecolor='lightgray', alpha=0.5))
    plt.gca().add_patch(plt.Rectangle((0.5, 0.5), 0.5, 0.5, facecolor='lightgray', alpha=0.5))

    # Create a legend with the specified labels or sequence numbers and corresponding colors
    if legend_labels:
        plt.legend(labels=legend_labels, loc='upper right')
    else:
        plt.legend(loc='upper right')

    # Show the plot
    plt.show()