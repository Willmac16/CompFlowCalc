def isentropic(gamma, mach):
  temp_ratio = 1 + (gamma - 1) / 2 * mach ** 2

  pressure_ratio = temp_ratio ** (gamma / (gamma - 1))
  density_ratio = temp_ratio ** (1 / (gamma - 1))
  return pressure_ratio, temp_ratio, density_ratio

def normal_shock(gamma, mach):
  pressure_ratio = (2 * gamma * mach ** 2 - (gamma - 1)) / (gamma + 1)
  temp_ratio = (2 * gamma * mach ** 2 - (gamma - 1)) * ((gamma - 1) * mach ** 2 + 2) / ((gamma + 1) ** 2 * mach ** 2)
  density_ratio = (gamma + 1) * mach ** 2 / ((gamma - 1) * mach ** 2 + 2)

  stag_ratio = temp_ratio ** (-gamma / (gamma - 1)) * pressure_ratio

  exit_stag_ratio = stag_ratio * (1 + (gamma - 1) / 2 * mach ** 2) ** (gamma / (gamma - 1))

  return pressure_ratio, temp_ratio, density_ratio, stag_ratio, exit_stag_ratio
