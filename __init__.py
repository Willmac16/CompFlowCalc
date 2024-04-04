import math
from scipy.optimize import root_scalar
import numpy as np

class IsentropicRatio:
    def __init__(self, gamma, mach):
        self.temp_ratio = 1 + (gamma - 1) / 2 * mach ** 2

        self.pressure_ratio = self.temp_ratio ** (gamma / (gamma - 1))
        self.density_ratio = self.temp_ratio ** (1 / (gamma - 1))
        self.mu = math.asin(1 / mach)

def isentropic(gamma, mach):
    return IsentropicRatio(gamma, mach)

class PressureInverseIsentropicRatio:
    def __init__(self, gamma, pressure_ratio):
        self.temp_ratio = pressure_ratio ** ((gamma - 1) / gamma)
        self.mach = (2 * (self.temp_ratio - 1) / (gamma - 1)) ** 0.5

def pressure_inverse_isentropic(gamma, pressure_ratio):
    return PressureInverseIsentropicRatio(gamma, pressure_ratio)

class NormalShockRatio:
    def __init__(self, gamma, mach):
        self.pressure_ratio = (2 * gamma * mach ** 2 - (gamma - 1)) / (gamma + 1)
        self.temp_ratio = (2 * gamma * mach ** 2 - (gamma - 1)) * ((gamma - 1) * mach ** 2 + 2) / ((gamma + 1) ** 2 * mach ** 2)
        self.density_ratio = (gamma + 1) * mach ** 2 / ((gamma - 1) * mach ** 2 + 2)

        self.stag_ratio = self.temp_ratio ** (-gamma / (gamma - 1)) * self.pressure_ratio

        stag_one = IsentropicRatio(gamma, mach).pressure_ratio
        self.exit_stag_ratio = self.stag_ratio * stag_one
        self.mach_two = ((1 + (gamma - 1) / 2 * mach ** 2) / (gamma * mach ** 2 - (gamma - 1) / 2)) ** 0.5

def normal_shock(gamma, mach):
    return NormalShockRatio(gamma, mach)


class ExpansionWaveRatio:
    def __init__(self, gamma, mach):
        self.velocity_ratio = 1 / (1 + (gamma - 1) / 2 * mach)
        self.pressure_ratio = 1 / (1 + (gamma - 1) / 2 * mach) ** (2 * gamma / (gamma - 1))
        self.temp_ratio = 1 / (1 + (gamma - 1) / 2 * mach) ** 2
        self.density_ratio = 1 / (1 + (gamma - 1) / 2 * mach) ** (2 / (gamma - 1))

def expansion_wave(gamma, mach):
    return ExpansionWaveRatio(gamma, mach)

def tbm_implicit(gamma, theta, beta, mach):
    lhs = math.tan(theta)
    rhs = 2 / math.tan(beta) * (mach ** 2 * math.sin(beta) ** 2 - 1) / (mach ** 2 * (gamma + math.cos(2 * beta)) + 2)
    return rhs - lhs

class ThetaBetaMach:
    def __init__(self, gamma, mach=math.nan, beta=math.nan, theta=math.nan):
        if not math.isnan(beta) and not math.isnan(mach):
            self.mach = mach
            self.beta = beta
            self.theta = math.atan(2 / math.tan(beta) * (mach ** 2 * math.sin(beta) ** 2 - 1) / (mach ** 2 * (gamma + math.cos(2 * beta)) + 2))
        elif not math.isnan(beta) and not math.isnan(theta):
            self.beta = beta
            self.theta = theta

            self.mach = root_scalar(lambda mach:
                                    tbm_implicit(gamma, theta, beta, mach),
                                    bracket=[0.1, 1000]).root

        elif not math.isnan(mach) and not math.isnan(theta):
            self.mach = mach
            self.theta = theta

            self.beta = root_scalar(lambda beta:
                                    tbm_implicit(gamma, theta, beta, mach),
                                    x0=7/9 * math.pi / 4,
                                    bracket=[1e-6, math.pi / 4],
                                    method='brentq').root

class ObliqueShockResult:
    def __init__(self, gamma, mach_one, beta):
        self.theta = math.atan(2 / math.tan(beta) * (mach_one ** 2 * math.sin(beta) ** 2 - 1) / (mach_one ** 2 * (gamma + math.cos(2 * beta)) + 2))
        mach_one_normal = mach_one * math.sin(beta)

        self.shock_ratio: NormalShockRatio = NormalShockRatio(gamma, mach_one_normal)
        mach_two_normal = self.shock_ratio.mach_two
        self.mach_two = mach_two_normal / math.sin(self.theta - beta)


class PrandtlMeyer:
    def __init__(self, gamma, mach):
        self.nu = math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (mach ** 2 - 1))) - math.atan(math.sqrt(mach ** 2 - 1))

if __name__ == "__main__":
    GAMMA = 1.4

    print("Theta @ M2, Beta = 50˚:", math.degrees(ThetaBetaMach(GAMMA, mach=2, beta=math.radians(50)).theta))

    assert abs(pressure_inverse_isentropic(GAMMA, isentropic(GAMMA, 2).pressure_ratio).mach - 2.0) < 1e-6

    print(f"Pressure ratio of a Mach 2 Normal Shock {normal_shock(GAMMA, 2).pressure_ratio:.5f}")
