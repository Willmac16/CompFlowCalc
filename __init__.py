class IsentropicRatio:
    def __init__(self, pressure_ratio, temp_ratio, density_ratio):
        self.pressure_ratio = pressure_ratio
        self.temp_ratio = temp_ratio
        self.density_ratio = density_ratio

def isentropic(gamma, mach):
    temp_ratio = 1 + (gamma - 1) / 2 * mach ** 2

    pressure_ratio = temp_ratio ** (gamma / (gamma - 1))
    density_ratio = temp_ratio ** (1 / (gamma - 1))
    return IsentropicRatio(pressure_ratio, temp_ratio, density_ratio)

class PressureInverseIsentropicRatio:
    def __init__(self, temp_ratio, mach):
        self.temp_ratio = temp_ratio
        self.mach = mach

def pressure_inverse_isentropic(gamma, pressure_ratio):
    temp_ratio = pressure_ratio ** ((gamma - 1) / gamma)

    mach = (2 * (temp_ratio - 1) / (gamma - 1)) ** 0.5
    return PressureInverseIsentropicRatio(temp_ratio, mach)

class NormalShockRatio:
    def __init__(self, pressure_ratio, temp_ratio, density_ratio, stag_ratio, exit_stag_ratio):
        self.pressure_ratio = pressure_ratio
        self.temp_ratio = temp_ratio
        self.density_ratio = density_ratio
        self.stag_ratio = stag_ratio
        self.exit_stag_ratio = exit_stag_ratio

def normal_shock(gamma, mach):
    pressure_ratio = (2 * gamma * mach ** 2 - (gamma - 1)) / (gamma + 1)
    temp_ratio = (2 * gamma * mach ** 2 - (gamma - 1)) * ((gamma - 1) * mach ** 2 + 2) / ((gamma + 1) ** 2 * mach ** 2)
    density_ratio = (gamma + 1) * mach ** 2 / ((gamma - 1) * mach ** 2 + 2)

    stag_ratio = temp_ratio ** (-gamma / (gamma - 1)) * pressure_ratio

    stag_one = isentropic(gamma, mach).pressure_ratio
    exit_stag_ratio = stag_ratio * stag_one

    return NormalShockRatio(pressure_ratio, temp_ratio, density_ratio, stag_ratio, exit_stag_ratio)


class ExpansionWaveRatio:
    def __init__(self, velocity_ratio, pressure_ratio, temp_ratio, density_ratio):
        self.velocity_ratio = velocity_ratio
        self.pressure_ratio = pressure_ratio
        self.temp_ratio = temp_ratio
        self.density_ratio = density_ratio

def expansion_wave(gamma, mach):
    velocity_ratio = 1 / (1 + (gamma - 1) / 2 * mach)
    pressure_ratio = 1 / (1 + (gamma - 1) / 2 * mach) ** (2 * gamma / (gamma - 1))
    temp_ratio = 1 / (1 + (gamma - 1) / 2 * mach) ** 2
    density_ratio = 1 / (1 + (gamma - 1) / 2 * mach) ** (2 / (gamma - 1))

    return ExpansionWaveRatio(velocity_ratio, pressure_ratio, temp_ratio, density_ratio)

if __name__ == "__main__":
    GAMMA = 1.4

    assert(abs(pressure_inverse_isentropic(GAMMA, isentropic(GAMMA, 2).pressure_ratio).mach - 2.0) < 1e-6)

    print("Pressure ratio of a Mach 2 Normal Shock {:.5f}".format(normal_shock(GAMMA, 2).pressure_ratio))
