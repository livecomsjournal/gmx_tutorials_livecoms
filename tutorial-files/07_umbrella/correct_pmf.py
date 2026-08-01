import numpy as np
import os

def apply_pmf_correction(input_file, output_file, temperature, bulk_cutoff):
	"""
	Applies spherical volume entropy correction to a GROMACS PMF profile.
	
	Parameters:
	input_file (str): Path to input profile.xvg
	output_file (str): Path for the corrected output xvg file
	temperature (float): Simulation temperature in Kelvin
	bulk_cutoff (float): Distance (r) above which the system is considered 'bulk'
	"""
	kB = 0.00831451  # kJ/(mol*K)
	factor = 2 * kB * temperature
	
	r_vals = []
	w_raw = []
	headers = []
	
	# Read the GROMACS xvg file
	with open(input_file, 'r') as f:
		for line in f:
			if line.startswith(('@', '#')):
				headers.append(line)
			else:
				parts = line.split()
				if len(parts) >= 2:
					r_vals.append(float(parts[0]))
					w_raw.append(float(parts[1]))
					
	r_vals = np.array(r_vals)
	w_raw = np.array(w_raw)
	
	# Avoid log(0) error if profile starts exactly at 0
	with np.errstate(divide='ignore', invalid='ignore'):
		correction = factor * np.log(r_vals)
		# Handle r=0 case gracefully if it exists
		correction[r_vals == 0] = 0.0
		
	w_corrected = w_raw + correction
	
	# Align bulk plateau to 0 kJ/mol
	bulk_indices = np.where(r_vals >= bulk_cutoff)[0]
	if len(bulk_indices) == 0:
		print(f"Warning: No points found above bulk cutoff {bulk_cutoff}. Using last 10 points.")
		bulk_mean = np.mean(w_corrected[-10:])
	else:
		bulk_mean = np.mean(w_corrected[bulk_indices])
		
	w_final = w_corrected - bulk_mean
	
	# Write corrected profile
	with open(output_file, 'w') as f:
		for header in headers:
			f.write(header)
		for r, w in zip(r_vals, w_final):
			f.write(f"{r:12.7f} {w:12.7f}\n")
			
	print(f"Successfully generated {output_file}")
	print(f"Shifted bulk baseline energy by subtracting: {bulk_mean:.4f} kJ/mol")

# --- USER CONFIGURATION ---
input_xvg = "profile.xvg"		 # Your raw PMF file from gmx wham
output_xvg = "pmf_corrected.xvg" # Output file name
temp_kelvin = 298.0				 # Simulation temperature
bulk_min_r = 2.0				 # Distance (nm) where interaction becomes zero

apply_pmf_correction(input_xvg, output_xvg, temp_kelvin, bulk_min_r)
