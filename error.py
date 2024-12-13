import re
import numpy as np

def parse_gaussian_data_from_file(file_path):
    """
    Reads the Gaussian output data from a file and extracts relevant parameters for error calculation.
    :param file_path: Path to the .txt file containing the Gaussian output data.
    :return: List of dictionaries with 'mu', 'sigma', and 'max_velocity' for each measurement point.
    """
    with open(file_path, 'r') as file:
        data = file.read()

    measurements = []
    # Split the data into blocks for each measurement
    data_blocks = data.split("# BEGINHEADER")
    for block in data_blocks:
        if "Peaks Estimates" in block:
            # Extract maximum probable velocity
            max_velocity_match = re.search(r"Maximum Probable Velocity = ([\d\.\-]+)", block)
            max_velocity = float(max_velocity_match.group(1)) if max_velocity_match else None

            # Extract peak parameters
            peaks = re.findall(
                r"Estimated centre for peak \d+ = ([\d\.\-]+)\s+Estimated std for peak \d+ = ([\d\.\-]+)",
                block
            )

            # Use the highest velocity peak (largest center value)
            if peaks:
                peaks = [(float(mu), float(sigma)) for mu, sigma in peaks]
                highest_peak = max(peaks, key=lambda x: x[0])  # Largest center value
                mu, sigma = highest_peak

                measurements.append({
                    'mu': mu,
                    'sigma': sigma,
                    'max_velocity': max_velocity
                })
    return measurements

def calculate_vr_errors(measurements, sigma_mu=0.5):
    """
    Calculates the errors for V_r,max based on the extracted Gaussian parameters.
    :param measurements: List of dictionaries with 'mu', 'sigma', and 'max_velocity'.
    :param sigma_mu: Assumed uncertainty in the peak center (default: 0.5 km/s).
    :return: List of errors for V_r,max.
    """
    errors = []
    for m in measurements:
        mu = m['mu']
        sigma = m['sigma']
        # Estimate sigma_sigma
        sigma_sigma = sigma / np.sqrt(256)  # Assuming 256 frequency channels
        # Propagate errors
        vr_error = np.sqrt(sigma_mu**2 + (3 * sigma_sigma)**2)
        errors.append(vr_error)
    return errors

# Example usage
file_path = "outputdata2.txt"  # Replace with the actual file path
measurements = parse_gaussian_data_from_file(file_path)
vr_errors = calculate_vr_errors(measurements)

# Print results
for i, (m, error) in enumerate(zip(measurements, vr_errors)):
    print(f"Measurement {i + 1}:")
    print(f"  Maximum Velocity: {m['max_velocity']} km/s")
    print(f"  Center: {m['mu']} km/s")
    print(f"  Std Dev: {m['sigma']} km/s")
    print(f"  Error in V_r,max: ±{error:.3f} km/s\n")
