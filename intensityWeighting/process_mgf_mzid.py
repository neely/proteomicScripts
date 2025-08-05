"""
Script: process_mgf_mzid.py

Description:
This script processes mass spectrometry data by integrating information from MGF, optional target CSV, and mzIdentML files.
It extracts scan data from MGF, performs optional m/z target matching, and incorporates identification
threshold status from mzIdentML, finally merging all relevant data into a single CSV output.

Input Files:
1.  <input.mgf>: The Mascot Generic Format (MGF) file containing MS/MS spectral data.
2.  <input.mzid>: The mzIdentML file, which provides spectrum identification results, including passThreshold status.
3.  <targets.csv> (Optional): A CSV file specifying target m/z values, their tolerances, and tolerance types (e.g., ppm or Da).

Command-Line Usage:
Run the script from your terminal using one of the following formats:

    python process_mgf_mzid.py <input.mgf> <input.mzid> <output.csv>
    (Processes MGF and mzIdentML files without target matching.)

    OR

    python process_mgf_mzid.py <input.mgf> <targets.csv> <input.mzid> <output.csv>
    (Processes MGF, targets.csv, and mzIdentML files, including target matching.)

Arguments:
- <input.mgf>: Path to the MGF data file.
- <targets.csv>: Path to the optional mass target CSV file.
- <input.mzid>: Path to the mzIdentML search results file.
- <output.csv>: Desired name/path for the resulting output CSV file.

Output File (<output.csv>):
The generated CSV file will contain the following columns for each processed scan:
- `scan_number`: Unique identifier for the MS2 scan from the MGF file.
- `RT`: Retention time (in seconds) of the scan.
- `pepmass`: Precursor m/z value of the scan.
- `intensity`: MS1 intensity of the precursor.
- `summed_ms2_intensity`: The sum of all MS2 ion intensities within the scan.
- `[target m/z]`: Columns (e.g., "123.456") indicating "yes" if an MS2 ion matched the specific target m/z, "no" otherwise. These columns are only present if a `<targets.csv>` file was provided.
- `pass_threshold`: Indicates "yes" if the scan passed the identification threshold in the mzIdentML file, "no" otherwise.

Author: B. Neely; 30 Jan 2025
"""


import csv
import sys
import xml.etree.ElementTree as ET


def parse_mass_targets(mass_targets_file):
    """Parses the target m/z CSV file and returns a list of target values with tolerances."""
    targets = []
    try:
        with open(mass_targets_file, 'r', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile)
            # Ensure required headers are present
            required_headers = ["m/z", "tolerance", "tolerance_type"]
            if not all(header in reader.fieldnames for header in required_headers):
                print(f"Error: Target CSV must contain all of these columns: {', '.join(required_headers)}")
                sys.exit(1)

            for row in reader:
                targets.append({
                    "mz": float(row["m/z"]),
                    "tolerance": float(row["tolerance"]),
                    "tolerance_type": row["tolerance_type"].strip().lower()
                })
        print(f"Successfully parsed {len(targets)} mass targets from {mass_targets_file}.")
    except FileNotFoundError:
        print(f"Error: Target file '{mass_targets_file}' not found.")
        sys.exit(1)
    except ValueError as e:
        print(f"Error parsing numerical data in '{mass_targets_file}': {e}. Ensure 'm/z' and 'tolerance' columns contain valid numbers.")
        sys.exit(1)
    except KeyError as e:
        print(f"Error: Missing expected column in '{mass_targets_file}': {e}. Please check your CSV headers.")
        sys.exit(1)
    return targets


def is_within_tolerance(mz_value, target, tolerance, tol_type):
    """Checks if a given m/z value is within the tolerance range of a target mass."""
    if tol_type == "ppm":
        tolerance_da = (tolerance / 1e6) * target  # Convert ppm to Da range
    else: # Assume Da if not ppm
        tolerance_da = tolerance
    return abs(mz_value - target) <= tolerance_da


def parse_mgf(mgf_file, targets):
    """Parses an MGF file and extracts scan data along with m/z target matching and MS2 summed intensity."""
    precursor_data = []
    print(f"Parsing MGF file: {mgf_file}...")
    try:
        with open(mgf_file, 'r') as f:
            inside_ion_block = False
            scan_number, rt, precursor_mz, intensity = None, None, None, None
            mz_values = []
            summed_ms2_intensity = 0  # Track summed MS2 intensity

            for line in f:
                line = line.strip()

                if line == "BEGIN IONS":
                    inside_ion_block = True
                    scan_number, rt, precursor_mz, intensity = None, None, None, None
                    mz_values = []
                    summed_ms2_intensity = 0  # Reset summed intensity for new ion block

                elif inside_ion_block:
                    if line.startswith("TITLE="):
                        if "scan=" in line:
                            # Extract scan number, handling various formats
                            try:
                                scan_number = line.split("scan=")[1].split('&quot;')[0].strip().replace('"', '').split(' ')[0]
                            except IndexError:
                                scan_number = "Unknown" # Fallback if scan format is unexpected
                    elif line.startswith("RTINSECONDS="):
                        try:
                            rt = float(line.split("=")[-1].strip())
                        except ValueError:
                            rt = None # Handle cases where RT is not a valid float
                    elif line.startswith("PEPMASS="):
                        pepmass_values = line.split("=")[-1].strip().split()
                        try:
                            precursor_mz = float(pepmass_values[0])
                            intensity = float(pepmass_values[1]) if len(pepmass_values) > 1 else 0.0
                        except ValueError:
                            precursor_mz = None
                            intensity = 0.0 # Handle cases where pepmass/intensity are not valid floats
                    elif line == "END IONS":
                        inside_ion_block = False
                        if scan_number and precursor_mz is not None:
                            # Initialize target matches as "no"
                            target_hits = {str(target["mz"]): "no" for target in targets}

                            # Check if any m/z values from MS2 spectrum match the targets
                            if targets: # Only perform target matching if targets were loaded
                                for mz, _ in mz_values:
                                    for target in targets:
                                        if is_within_tolerance(mz, target["mz"], target["tolerance"], target["tolerance_type"]):
                                            target_hits[str(target["mz"])] = "yes"

                            precursor_data.append({
                                "scan_number": scan_number,
                                "RT": rt,
                                "pepmass": precursor_mz,
                                "intensity": intensity,  # MS1 Intensity
                                "summed_ms2_intensity": summed_ms2_intensity,  # Sum of MS2 ion intensities
                                **target_hits  # Columns for each target mass (if targets provided)
                            })
                    else: # Process MS2 ion data (m/z intensity pairs)
                        parts = line.split()
                        if len(parts) == 2:
                            try:
                                mz, intensity_value = float(parts[0]), float(parts[1])
                                mz_values.append((mz, intensity_value))
                                summed_ms2_intensity += intensity_value
                            except ValueError:
                                pass  # Ignore lines that aren't valid m/z-intensity pairs
        print(f"Finished parsing MGF file. Found {len(precursor_data)} scans.")
    except FileNotFoundError:
        print(f"Error: MGF file '{mgf_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while parsing MGF file '{mgf_file}': {e}")
        sys.exit(1)
    return precursor_data


def parse_mzid_with_scans(mzid_file):
    """Parses an mzIdentML file and extracts scan numbers with passThreshold status."""
    mzid_data = {}
    # Define namespaces to correctly parse XML
    namespaces = {
        'ns': 'http://psidev.info/psi/pi/mzIdentML/1.1',
        'xsi': 'http://www.w3.org/2001/XMLSchema-instance' # Often needed for schemaLocation
    }

    print(f"Parsing mzIdentML file: {mzid_file}...")
    try:
        tree = ET.parse(mzid_file)
        root = tree.getroot()

        # Iterate through SpectrumIdentificationResult elements
        for sir in root.findall('.//ns:SpectrumIdentificationResult', namespaces):
            scan_number = None
            # Extract spectrum title and parse scan number
            spectrum_title_elem = sir.find('./ns:cvParam[@name="spectrum title"]', namespaces)
            if spectrum_title_elem is not None:
                value = spectrum_title_elem.attrib.get('value', '')
                if 'scan=' in value:
                    try:
                        # Improved parsing for scan number from title
                        scan_number = value.split('scan=')[1].split('&quot;')[0].strip().replace('"', '').split(' ')[0]
                    except IndexError:
                        scan_number = "Unknown" # Fallback

            # Determine passThreshold status
            # A SpectrumIdentificationItem with passThreshold="true" indicates the scan passed
            pass_threshold = sir.find('./ns:SpectrumIdentificationItem[@passThreshold="true"]', namespaces) is not None

            if scan_number is not None:
                mzid_data[scan_number] = pass_threshold
        print(f"Finished parsing mzIdentML file. Found {len(mzid_data)} identified scans.")
    except FileNotFoundError:
        print(f"Error: mzIdentML file '{mzid_file}' not found.")
        sys.exit(1)
    except ET.ParseError as e:
        print(f"Error parsing mzIdentML file '{mzid_file}': {e}. Ensure it's a valid XML file.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while parsing mzIdentML file '{mzid_file}': {e}")
        sys.exit(1)
    return mzid_data


def merge_data(mgf_data, mzid_data):
    """Merges MGF and mzID data by scan number, adding pass_threshold info."""
    print("Merging MGF and mzIdentML data...")
    for entry in mgf_data:
        scan_number = entry["scan_number"]
        # Default to "no" if scan not found in mzIdentML data
        entry["pass_threshold"] = "yes" if mzid_data.get(scan_number, False) else "no"
    print("Data merge complete.")
    return mgf_data


def export_to_csv(data, output_file, targets):
    """Exports the processed data to a CSV file, including summed MS2 intensity and optional target columns."""
    print(f"Exporting data to {output_file}...")
    if not data:
        print("No data to export.")
        return

    # Dynamically build fieldnames based on whether targets were provided
    base_fieldnames = ["scan_number", "RT", "pepmass", "intensity", "summed_ms2_intensity"]
    target_fieldnames = [str(t["mz"]) for t in targets] if targets else []
    final_fieldnames = base_fieldnames + target_fieldnames + ["pass_threshold"]

    try:
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=final_fieldnames)
            writer.writeheader()
            for row in data:
                # Ensure all fields in final_fieldnames are present in the row, add empty string if not
                # This handles cases where target_hits might not be present if targets were not used
                cleaned_row = {key: row.get(key, '') for key in final_fieldnames}
                writer.writerow(cleaned_row)
        print(f"Data successfully exported to {output_file}.")
    except IOError as e:
        print(f"Error writing to output file '{output_file}': {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during CSV export: {e}")
        sys.exit(1)


def main():
    # Adjust argument check for optional targets.csv
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python process_mgf_mzid.py <input.mgf> <input.mzid> <output.csv>")
        print("OR:    python process_mgf_mzid.py <input.mgf> <targets.csv> <input.mzid> <output.csv>")
        sys.exit(1)

    # Determine which arguments are which based on count
    if len(sys.argv) == 4:
        mgf_file = sys.argv[1]
        mzid_file = sys.argv[2]
        output_csv = sys.argv[3]
        targets = [] # No targets.csv provided, so targets list is empty
        print("Running without mass targets file.")
    else: # len(sys.argv) == 5
        mgf_file = sys.argv[1]
        targets_file = sys.argv[2]
        mzid_file = sys.argv[3]
        output_csv = sys.argv[4]
        targets = parse_mass_targets(targets_file) # Parse targets.csv
        print(f"Running with mass targets file: {targets_file}.")


    mgf_data = parse_mgf(mgf_file, targets)
    mzid_data = parse_mzid_with_scans(mzid_file)
    merged_data = merge_data(mgf_data, mzid_data)

    export_to_csv(merged_data, output_csv, targets)


if __name__ == "__main__":
    main()
