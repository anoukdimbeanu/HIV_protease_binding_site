import RNA

# Define your RNA sequence
sequence = "GCGCUUCGCCG"

# Create a fold compound object
fc = RNA.fold_compound(sequence)

# Get minimum free energy structure
(ss, mfe) = fc.mfe()

print(f"RNA Sequence: {sequence}")
print(f"Secondary Structure: {ss}")
print(f"Minimum Free Energy: {mfe} kcal/mol")

# Calculate base pair probabilities
fc.pf()  # partition function needed before accessing probabilities

# Print some pair probabilities
print("\nBase Pair Probabilities (i, j, p):")
length = len(sequence)
for i in range(1, length + 1):
    for j in range(i + 1, length + 1):
        p = fc.bpp(i, j)
        if p > 0.01:
            print(f"{i}-{j}: {p:.4f}")

# Optional: Dot plot with matplotlib
try:
    import matplotlib.pyplot as plt
    import numpy as np

    # Create dot plot matrix
    dot_plot = np.zeros((length, length))
    for i in range(1, length + 1):
        for j in range(i + 1, length + 1):
            p = fc.bpp(i, j)
            if p > 0.01:
                dot_plot[i-1][j-1] = p

    # Plot
    plt.imshow(dot_plot, cmap="hot", interpolation="nearest")
    plt.title("Base Pair Probability Dot Plot")
    plt.xlabel("Position")
    plt.ylabel("Position")
    plt.colorbar(label="Pairing Probability")
    plt.show()

except ImportError:
    print("\nInstall matplotlib and numpy to plot dot plots:")
    print("pip install matplotlib numpy")
