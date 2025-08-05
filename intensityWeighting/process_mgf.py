"""
Script: process_mgf.py

Description:
This script processes mass spectrometry data from an MGF (Mascot Generic Format) file.
It extracts key scan information, calculates the summed intensity of MS2 ions, and
optionally checks for matches against a list of predefined target m/z values.
The processed data is then exported to a CSV file.

Input Files:
1.  <input.mgf>: The MGF file containing MS/MS spectral data.
2.  <mass-targets.csv> (Optional): A CSV file specifying target m/z values, their tolerances,
    and tolerance types (e.g., ppm or Da).

Command-Line Usage:
Run the script from your terminal using one of the following formats:

    python process_mgf.py <input.mgf> <output.csv>
    (Processes MGF file without target m/z matching.)

    OR

    python process_mgf.py <input.mgf> <output.csv> <mass-targets.csv>
    (Processes MGF file and performs target m/z matching.)

Arguments:
- <input.mgf>: Path to the MGF data file.
- <output.csv>: Desired name/path for the resulting output CSV file.
- <mass-targets.csv>: Path to the optional mass target CSV file.

Output File (<output.csv>):
The generated CSV file will contain the following columns for each processed scan:
- `scan_number`: Unique identifier for the MS2 scan from the MGF file.
- `RT`: Retention time (in seconds) of the scan.
- `pepmass`: Precursor m/z value of the scan.
- `intensity`: MS1 intensity of the precursor.
- `MS2_sum`: The sum of all MS2 ion intensities within the scan.
- `[target m/z]`: Columns (e.g., "123.456") indicating "yes" if an MS2 ion matched the specific target m/z, "no" otherwise. These columns are only present if a `<mass-targets.csv>` file was provided.

Author: B. Neely; 30 Jan 2025
"""

import csv
import sys

def parse_mass_targets(target_file):
    """Parse the mass targets CSV file, handling potential BOM issues."""
    with open(target_file, 'r', encoding='utf-8-sig') as f:  # <-- This automatically removes BOM
        reader = csv.DictReader(f)
        print("Headers found:", reader.fieldnames)  # Debugging line

        targets = []
        for row in reader:
            print("Row data:", row)  # Debugging line
            targets.append({
                "mz": float(row["m/z"].strip()),  # Now it should work!
                "tolerance": float(row["tolerance"].strip()),
                "tolerance_type": row["tolerance_type"].strip().lower()
            })
    return targets

def is_within_tolerance(mz_value, target, tolerance, tol_type):
    """Checks if a given m/z value is within the tolerance range of a target mass."""
    if tol_type == "ppm":
        tolerance = (tolerance / 1e6) * target  # Convert ppm to Da range
    return abs(mz_value - target) <= tolerance

def parse_mgf(mgf_file, targets):
    """Parses an MGF file and extracts scan number, RT, precursor intensity, MS2 sum, and target matches."""
    precursor_data = []
    
    print("Parsing MGF file...")
    with open(mgf_file, 'r') as f:
        inside_ion_block = False
        scan_number = None
        precursor_mz = None
        intensity = None
        rt = None
        ms2_sum = 0
        mz_values = []

        for line in f:
            line = line.strip()

            # Identify start of an ion block
            if line == "BEGIN IONS":
                inside_ion_block = True
                scan_number = None
                precursor_mz = None
                intensity = None
                rt = None
                ms2_sum = 0
                mz_values = []

            elif inside_ion_block:
                if line.startswith("TITLE="):
                    if "scan=" in line:
                        scan_number = line.split("scan=")[-1].strip().split()[0].strip('"')
                elif line.startswith("RTINSECONDS="):
                    rt = float(line.split("=")[-1].strip())
                elif line.startswith("PEPMASS="):
                    pepmass_values = line.split("=")[-1].strip().split()
                    precursor_mz = float(pepmass_values[0])
                    intensity = float(pepmass_values[1]) if len(pepmass_values) > 1 else 0
                elif line == "END IONS":
                    # Check target masses
                    target_hits = {str(target["mz"]): "no" for target in targets}
                    for mz, _ in mz_values:
                        for target in targets:
                            if is_within_tolerance(mz, target["mz"], target["tolerance"], target["tolerance_type"]):
                                target_hits[str(target["mz"])] = "yes"

                    # Store the parsed data
                    if scan_number and precursor_mz is not None:
                        precursor_data.append({
                            'scan_number': scan_number,
                            'RT': rt,
                            'pepmass': precursor_mz,
                            'intensity': intensity,
                            'MS2_sum': ms2_sum,
                            **target_hits
                        })
                    inside_ion_block = False

                else:
                    # Process MS2 ion data
                    parts = line.split()
                    if len(parts) == 2:
                        try:
                            mz, intensity_val = float(parts[0]), float(parts[1])
                            mz_values.append((mz, intensity_val))
                            ms2_sum += intensity_val
                        except ValueError:
                            pass  # Ignore lines that aren't m/z-intensity pairs

    return precursor_data

def export_to_csv(data, output_file):
    """Exports parsed data to a CSV file."""
    print(f"Exporting data to {output_file}...")
    if not data:
        print("No data to export.")
        return

    fieldnames = list(data[0].keys())  # Ensure all columns are written
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

    print(f"Data successfully exported to {output_file}.")

def main():
    if len(sys.argv) < 3:
        print("Usage: python process_mgf.py <input.mgf> <output.csv> [mass-targets.csv]")
        sys.exit(1)

    mgf_file = sys.argv[1]
    output_csv = sys.argv[2]
    mass_targets_file = sys.argv[3] if len(sys.argv) > 3 else None

    targets = parse_mass_targets(mass_targets_file) if mass_targets_file else []
    precursor_data = parse_mgf(mgf_file, targets)
    export_to_csv(precursor_data, output_csv)

if __name__ == "__main__":
    main()
