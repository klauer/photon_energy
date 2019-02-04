from collections import namedtuple
import numpy as np

import matplotlib.pyplot as plt


def calculate_k(k0, a, b, period, gap):
    '''
    Parameters
    ----------
    k0 : float
        Tuning parameter
    a : float
        Tuning parameter
    b : float
        Tuning parameter
    period : float
        Undulator period, in mm
    gap : float, np.array
        Undulator gap, in mm
    '''
    return k0 * np.exp(a * (gap / period) + b * (gap / period) ** 2)


def calculate_photon_energy(electron_energy, period, k):
    '''
    Parameters
    ----------
    electron_energy : float
        Electron energy in GeV
    period : float
        Undulator period, in mm
    k : float
        Undulator K parameter (strength), unitless
    '''
    m_e = 0.0005109989461  # GeV [electron rest mass/energy]
    h = 6.62607004e-34  # Js [Planck's constant]
    e = 1.6021766208e-19  # C [electron charge]
    c = 299792458  # m/s, speed of light
    return ((2. * (electron_energy / m_e) ** 2 * h * c) /
            (e * period * 1e-3 * (1 + k ** 2 / 2.)))


UndulatorParameters = namedtuple('UndulatorParameters',
                                 'k0 a b period gaps electron_energies')

params = {
    'SXU': UndulatorParameters(k0=13.997,
                               a=-5.131,
                               b=1.878,
                               period=39.,
                               gaps=np.arange(7.2, 20, 0.1),
                               electron_energies=np.arange(3.6, 4.0, 0.001),
                               ),
    # hxu assuming copper line (SCRF line is slightly different, gap-wise)
    'HXU': UndulatorParameters(k0=9.471,
                               a=-5.131,
                               b=1.878,
                               period=26.,
                               gaps=np.arange(7.2, 19, 0.1),
                               electron_energies=np.arange(2.5, 15.0, 0.001),
                               ),
}


if __name__ == '__main__':
    for figure_no, line in enumerate(('HXU', 'SXU')):
        plt.figure(figure_no)
        k0, a, b, period, gaps, e_e = params[line]
        k = calculate_k(k0, a, b, period, gaps)

        for e_e0 in e_e[::len(e_e) // 10]:
            plt.figure(figure_no)
            e_r = calculate_photon_energy(e_e0, period, k)
            plt.plot(gaps, e_r, label='e_e={:.3f}GeV'.format(e_e0))
            plt.title(f'{line} Photon Energy vs Gap')
            plt.ylabel('Photon energy [eV]')
            plt.xlabel('Gap [mm]')
            plt.legend()

            plt.figure(2 + figure_no)
            plt.plot(k, e_r, label='e_e={:.3f}GeV'.format(e_e0))
            plt.title(f'{line} Photon Energy vs K')
            plt.ylabel('Photon energy [eV]')
            plt.xlabel('Undulator strength K')
            plt.legend()

    plt.show()
