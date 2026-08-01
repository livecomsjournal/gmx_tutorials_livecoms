import math

def compute_5_restraint_correction(r0, theta_A_deg, theta_B_deg, kr, k_thetaA, k_thetaB, k_phiA, k_phiB, T=298.15):
    """
    Computes the free energy correction for removing a 5-restraint setup 
    (1 distance, 2 angles, 2 dihedrals) using standard GROMACS units.
    
    Parameters:
    r0          : float -> Equilibrium distance (nm)
    theta_A_deg : float -> Equilibrium angle A (degrees)
    theta_B_deg : float -> Equilibrium angle B (degrees)
    kr          : float -> Distance force constant (kJ/mol/nm^2)
    k_thetaA    : float -> Angle A force constant (kJ/mol/rad^2)
    k_thetaB    : float -> Angle B force constant (kJ/mol/rad^2)
    k_phiA      : float -> Dihedral A force constant (kJ/mol/rad^2)
    k_phiB      : float -> Dihedral B force constant (kJ/mol/rad^2)
    T           : float -> Temperature (Kelvin), default is 298.15 K
    
    Returns:
    float: dG correction in kJ/mol
    """
    # 1. Physical Constants
    kB = 0.008314462618  # Boltzmann constant in kJ/(mol*K)
    V0 = 1.66054         # Standard state volume (1 M) in nm^3
    
    # 2. Convert Angles from Degrees to Radians
    theta_A = math.radians(theta_A_deg)
    theta_B = math.radians(theta_B_deg)
    
    # 3. Compute Phase-Space Volume Prefactor
    numerator = 2 * math.pi * V0
    denominator = (r0**2) * math.sin(theta_A) * math.sin(theta_B)
    prefactor = numerator / denominator
    
    # 4. Compute Harmonic Force Constant Matrix Component
    fc_product = kr * k_thetaA * k_thetaB * k_phiA * k_phiB
    denom_sqrt = (2 * math.pi * kB * T)**5
    sqrt_term = math.sqrt(fc_product / denom_sqrt)
    
    # 5. Calculate Final Free Energy Change
    inside_ln = prefactor * sqrt_term
    dG_correction = -kB * T * math.log(inside_ln)
    
    return dG_correction

# --- EXAMPLE SYSTEM WORKLOAD ---
# Define example parameters (adjust these to match your topology/mdp files)
params = {
    "r0": 0.53,           # 0.40 nm (4.0 Angstroms)
    "theta_A_deg": 95.57, # 90 degrees
    "theta_B_deg": 74.45, # 75 degrees
    "kr": 500.0,          # kJ/mol/nm^2
    "k_thetaA": 500.0,    # kJ/mol/rad^2
    "k_thetaB": 500.0,    # kJ/mol/rad^2
    "k_phiA": 500.0,      # kJ/mol/rad^2
    "k_phiB": 500.0,      # kJ/mol/rad^2
    "T": 298.15           # Kelvin
}

dg = compute_5_restraint_correction(**params)

print("--- 5-Restraint Correction Results ---")
print(f"Temperature  : {params['T']} K")
print(f"ΔG Correction: {dg:.4f} kJ/mol")
print(f"ΔG Correction: {dg / 4.184:.4f} kcal/mol")
