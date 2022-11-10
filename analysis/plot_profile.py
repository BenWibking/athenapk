import numpy as np
import matplotlib.pyplot as plt

code_length_cgs = 3.085677580962325e+24 # 1 Mpc
code_mass_cgs = 1.98841586e+47          # 1e14 Msun
code_time_cgs = 3.15576e+16             # 1 Gyr

if __name__ == '__main__':
    z, P, K, rho, n, ne, T, g, dP_dz = np.loadtxt("test_he_box.dat", unpack=True)
    
    code_density_cgs = code_mass_cgs / code_length_cgs**3

    plt.figure(figsize=(4,4), dpi=300)
    plt.plot(z, rho * code_density_cgs, label='density')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-3, 1)
    plt.xlabel(r'height $z$ (Mpc)')
    plt.ylabel(r'density $\rho$ (g cm$^{-3}$)')
    plt.tight_layout()
    plt.savefig("profile_rho.png")

    plt.figure(figsize=(4,4), dpi=300)
    plt.plot(z, T, label='temperature')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-3, 1)
    plt.ylim(1e6, 1e8)
    plt.xlabel(r'height $z$ (Mpc)')
    plt.ylabel(r'temperature (K)')
    plt.tight_layout()
    plt.savefig('profile_T.png')

    plt.figure(figsize=(4,4), dpi=300)
    plt.plot(z, K, label='entropy function $K$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-3, 1)
    plt.xlabel(r'height $z$ (Mpc)')
    plt.ylabel(r'entropy function $K$')
    plt.tight_layout()
    plt.savefig('profile_K.png')

    plt.figure(figsize=(4,4), dpi=300)
    plt.plot(z, g)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-3, 1)
    plt.xlabel(r'height $z$ (Mpc)')
    plt.ylabel(r'gravitational acceleration $|g|$')
    plt.tight_layout()
    plt.savefig('profile_g.png')

    plt.figure(figsize=(4,4), dpi=300)
    plt.plot(z, np.abs(dP_dz))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-3, 1)
    plt.xlabel(r'height $z$ (Mpc)')
    plt.ylabel(r'pressure gradient $dP/dz$')
    plt.tight_layout()
    plt.savefig('profile_dP_dz.png')
    