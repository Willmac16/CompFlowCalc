def isentropic(gamma, mach):
    temp_ratio = 1 + (gamma - 1) / 2 * mach ** 2

    pressure_ratio = temp_ratio ** (gamma / (gamma - 1))
    density_ratio = temp_ratio ** (1 / (gamma - 1))
    return pressure_ratio, temp_ratio, density_ratio

def press_inverse_isentropic(gamma, pressure_ratio):
    temp_ratio = pressure_ratio ** ((gamma - 1) / gamma)

    mach = (2 * (temp_ratio - 1) / (gamma - 1)) ** 0.5
    return temp_ratio, mach

def normal_shock(gamma, mach):
    pressure_ratio = (2 * gamma * mach ** 2 - (gamma - 1)) / (gamma + 1)
    temp_ratio = (2 * gamma * mach ** 2 - (gamma - 1)) * ((gamma - 1) * mach ** 2 + 2) / ((gamma + 1) ** 2 * mach ** 2)
    density_ratio = (gamma + 1) * mach ** 2 / ((gamma - 1) * mach ** 2 + 2)

    stag_ratio = temp_ratio ** (-gamma / (gamma - 1)) * pressure_ratio

    stag_one, _, _ = isentropic(gamma, mach)
    exit_stag_ratio = stag_ratio * stag_one

    return pressure_ratio, temp_ratio, density_ratio, stag_ratio, exit_stag_ratio
